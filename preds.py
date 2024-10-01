import os
import sys
import json
import time
import torch
import random
import argparse
import numpy as np
from torch import nn
from tqdm import tqdm
from math import ceil
from torch import FloatTensor as FT
from torch.autograd import Variable as V
from scunet import SCUNet
from utils import parse_map, write_map, split_map_into_overlapped_chunks, get_map_from_overlapped_chunks
import warnings
warnings.filterwarnings('ignore')

def pjoin(*args):
    return os.path.join(*args)

def softmax(x, axis=0):
    exp_x = np.exp(x - np.max(x, axis=axis, keepdims=True))
    return exp_x / np.sum(exp_x, axis=axis, keepdims=True)

def seed_torch(seed=42):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    # if you are using multi-GPU.
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.enabled = True
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"

seed_torch(42)

def get_args():
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", "-i", type=str, required=True, help="Input EM density map file")
    parser.add_argument("--output", "-o", type=str, required=True, help="Output directory of predicted")
    parser.add_argument("--contour", "-c", type=float, default=1e-6, help="Input contour level")
    parser.add_argument("--batchsize", "-b", type=int, default=40, help="Batchsize for prediction")
    parser.add_argument("--gpuid", "-g", type=str, help="Which GPU to use, '0' for #0", default='0')
    parser.add_argument("--models", "-m", type=str, help="Directory to deep learning models", default='./')
    parser.add_argument("--stride", "-s", type=int, help="Stride for splitting chunks", default=12)
    parser.add_argument("--usecpu", action='store_true', help="Run prediction on CPU")
    # Optional interp
    parser.add_argument("--pyinterp", action='store_true', help="Interpolation using pyinterp3d")
    args = parser.parse_args()
    return args


def load_model_and_run_inference_on_map(model_file, map_file, **kwargs):
    # Load map data
    apix = kwargs['apix']
    stride = kwargs['stride']
    box_size = kwargs['box_size']
    n_classes = kwargs['n_classes']
    use_gpu = kwargs['use_gpu']
    batch_size = kwargs['batch_size']
    gpu_id = kwargs['gpu_id']
    interp = kwargs['interp']

    n_gpus = 0
    if use_gpu:
        #os.environ["CUDA_VISIBLE_DEVICES"] = gpu_id
        if torch.cuda.is_available():
            n_gpus = torch.cuda.device_count()

    print("# Load map data from {}".format(map_file))
    map, origin, nxyz, voxel_size = parse_map(map_file, ignorestart=False, apix=apix, interp=interp)
    # Clip stride lower-bound for too large maps
    if np.min(nxyz) > 400:
        stride = max(stride, 16)
    print("# Map dimensions = {}".format(nxyz))
    print("# Using stride = {}".format(stride))
    chunks, ncx, ncy, ncz = split_map_into_overlapped_chunks(map, box_size, stride, dtype="float32", padding=0.0)
    n_chunks = len(chunks)
    print("# Split map into {} chunk(s)".format(n_chunks))

    # Normalize chunks and del non-positive chunks
    del_indices = []
    maximum = np.percentile(map[map > 0], 99.999)
    chunks_norm = np.zeros((n_chunks, box_size, box_size, box_size), dtype=np.float32)
    for i, chunk in enumerate(chunks):
        if chunk.max() <= 0.0:
            del_indices.append(i)
            continue
        chunks_norm[i] = chunk
    chunks_norm = chunks_norm.clip(min=0.0, max=maximum) / maximum * 100
    chunks = np.delete(chunks_norm, del_indices, axis=0)
    keep_indices = np.delete(np.arange(n_chunks, dtype=np.int32), del_indices, axis=0)
    n_chunks0 = len(chunks)
    print("# Get {} positive chunks".format(n_chunks0))
    X = V(FT(chunks), requires_grad=False).view(-1, 1, box_size, box_size, box_size) 
    del chunks
    data_num = len(X)
    n_batches = ceil(data_num/batch_size)

    # Load pretrained model stage1
    print("# Load trained weights from {}".format(model_file))
    if use_gpu:
        model = torch.load(model_file)
    else:
        model = torch.load(model_file, map_location=torch.device('cpu'))
    model_state_dict = model.module.state_dict()
    del model
    model = SCUNet(input_resolution=box_size, n_classes=n_classes)
    model.load_state_dict(model_state_dict)
    if use_gpu:
        torch.cuda.empty_cache()
        model = model.cuda()
        if n_gpus > 1:
            model = nn.DataParallel(model)
    model.eval()

    # Run inference
    chunks_pred0 = np.zeros((n_chunks0, n_classes, box_size, box_size, box_size), dtype="float32")
    print("# Start prediction", flush=True)
    with torch.no_grad():
        for i in tqdm(range(n_batches), ascii=True, ncols=50):
            #if i % 50 == 0:
            #    print("# Batch {}/{}".format(i+1, n_batches))

            X_batch = X[i * batch_size : (i + 1) * batch_size]
            if use_gpu:
                X_batch = X_batch.cuda()
            y_pred = model(X_batch) # [B, n_classes, L, L, L]
            y_pred0 = y_pred.cpu().detach().numpy()
            chunks_pred0[i * batch_size : (i + 1) * batch_size] = y_pred0
    chunks_pred = np.zeros((n_chunks, n_classes, box_size, box_size, box_size), dtype="float32")
    chunks_pred[keep_indices] = chunks_pred0
    del chunks_pred0
    map_pred0 = get_map_from_overlapped_chunks(chunks_pred, ncx, ncy, ncz, n_classes, box_size, stride, nxyz) # [n_classes, nx, ny, nz]
    return map_pred0, map, origin, nxyz, voxel_size

# Data params
data_params = {
    "apix": 1.0,
    "box_size": 48,
    "stride": 12,
    "interp": 'f90',
}

# Test params
test_params = {
    "use_gpu": True,
    "gpu_id": "0",
    "batch_size": 160
}

# Only inference stage1
def inference_stage1(dir_map, contour, dir_models, dir_out, data_params, test_params):

    print(f"# Select map contour at {contour:.6f}", flush=True)

    # Read hyper-params
    use_gpu = test_params['use_gpu']
    gpu_id = test_params['gpu_id']
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    if use_gpu:
        if gpu_id != 'auto':
            os.environ["CUDA_VISIBLE_DEVICES"] = gpu_id
        if torch.cuda.is_available():
            n_gpus = torch.cuda.device_count()
            print("# Running on {} GPU(s).".format(n_gpus))
            print("# Running on GPU {}.".format(gpu_id))
        else:
            print("CUDA not available.")
            sys.exit()
    else:
        n_gpus = 0
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        print("# Running on CPU.")

    batch_size = test_params['batch_size']

    apix = data_params['apix']
    box_size = data_params['box_size']
    assert box_size % 2 == 0
    stride = data_params['stride']

    # Run inference on stage1
    map_pred, map, origin, _, voxel_size = load_model_and_run_inference_on_map(
        model_file=os.path.join(dir_models, "stage1"),
        map_file=dir_map,

        apix=apix,    
        box_size=box_size,
        stride=stride,
        n_classes=3,
        use_gpu=use_gpu,
        batch_size=batch_size,
        gpu_id=gpu_id,
        interp=data_params['interp'],
    )

    # Write map of stage1
    map_pred = np.argmax(map_pred, axis=0)
    below_contour = np.where(map <= contour, 3, 0)
    map_pred = np.where(below_contour > map_pred, below_contour, map_pred)
    types = ["prot.mrc", "na.mrc", "bg.mrc"]
    # Only write na density map
    for i in [1]:
        mask = np.where(map_pred == i, 1, 0)
        out = mask * map
        dir_map_out = os.path.join(dir_out, types[i])
        write_map(dir_map_out, out.astype(np.float32), voxel_size, origin=origin)
        print("# Write map to {}".format(dir_map_out), flush=True)

    # Do a quick statistic on the ratio of NA/PROT/BG region
    n_prot = np.where(map_pred == 0, 1, 0).sum()
    n_na   = np.where(map_pred == 1, 1, 0).sum()
    n_bg   = np.where(map_pred == 2, 1, 0).sum()

    print("# Region report")
    print("# Among all voxels. prot ratio is {:.4f}".format( n_prot / (n_prot + n_na + n_bg + 1e-3) ))
    print("# Among all voxels. na   ratio is {:.4f}".format( n_na / (n_prot + n_na + n_bg + 1e-3) ))
    print("# Among all voxels. bg   ratio is {:.4f}".format( n_bg / (n_prot + n_na + n_bg + 1e-3) ))
    print("# Among prot and na voxels. prot ratio is {:.4f}".format( n_prot / (n_prot + n_na + 1e-3) ))
    print("# Among prot and na voxels. na   ratio is {:.4f}".format( n_na / (n_prot + n_na + 1e-3) ), flush=True)


# Only inference stage2
def inference_stage2(dir_map, contour, dir_models, dir_out, data_params, test_params):

    print(f"# Select map contour at {contour:.6f}", flush=True)

    # Read hyper-params
    use_gpu = test_params['use_gpu']
    gpu_id = test_params['gpu_id']
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    if use_gpu:
        if gpu_id != 'auto':
            os.environ["CUDA_VISIBLE_DEVICES"] = gpu_id
        if torch.cuda.is_available():
            n_gpus = torch.cuda.device_count()
            print("# Running on {} GPU(s).".format(n_gpus))
            print("# Running on GPU {}.".format(gpu_id))
        else:
            print("CUDA not available.")
            sys.exit()
    else:
        n_gpus = 0
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        print("# Running on CPU.")

    batch_size = test_params['batch_size']

    apix = data_params['apix']
    box_size = data_params['box_size']
    assert box_size % 2 == 0
    stride = data_params['stride']

    map_pred, map, origin, _, voxel_size = load_model_and_run_inference_on_map(
        model_file=os.path.join(dir_models, "stage2_prob"),
        map_file=dir_map,

        apix=apix,
        box_size=box_size,
        stride=stride,
        n_classes=3,
        use_gpu=use_gpu,
        batch_size=batch_size,
        gpu_id=gpu_id,
        interp=data_params['interp'],
    )
    mask = np.where(map <= contour, 0, 1).astype(np.int8)
    types = ["p.mrc", "c4.mrc", "n.mrc"]
    for i in range(3):
        out = mask * map_pred[i]
        dir_map_out = os.path.join(dir_out, types[i])
        write_map(dir_map_out, out.astype(np.float32), voxel_size, origin=origin)
        print("# Write map to {}".format(dir_map_out), flush=True)
    
    #map_merged = (map_pred[0] + map_pred[1]) / 2.
    #out = mask * map_merged
    #dir_map_out = os.path.join(dir_out, "merged.mrc")
    #write_map(dir_map_out, out.astype(np.float32), voxel_size, origin=origin)
    #print("# Write map to {}".format(dir_map_out), flush=True)


    map_pred, map, origin, _, voxel_size = load_model_and_run_inference_on_map(
        model_file=os.path.join(dir_models, "stage2_aa"),
        map_file=dir_map,

        apix=apix,
        box_size=box_size,
        stride=stride,
        n_classes=4,
        use_gpu=use_gpu,
        batch_size=batch_size,
        gpu_id=gpu_id,
        interp=data_params['interp'],
    )
    dir_map_out = os.path.join(dir_out, "prob.npz")
    prob_pred = softmax(map_pred, axis=0)
    np.savez(dir_map_out, prob=prob_pred.astype(np.float32))
    print("# Save NA probability npz to {}".format(dir_map_out))

    map_pred = np.argmax(map_pred, axis=0)
    below_contour = np.where(map <= contour, 4, 0)
    map_pred = np.where(below_contour > map_pred, below_contour, map_pred)
    dir_map_out = os.path.join(dir_out, "aa.mrc")
    out = map_pred
    write_map(dir_map_out, out.astype(np.float32), voxel_size, origin=origin)
    print("# Write map to {}".format(dir_map_out), flush=True)


if __name__ == "__main__":
    start = time.time()
    args = get_args()

    # Do deep learning prediction
    dir_map = args.input
    dir_out = args.output
    contour = args.contour
    dir_models = args.models
    contour = args.contour

    print("# Making directory {}".format(dir_out), flush=True)
    os.makedirs(dir_out, exist_ok=True)

    if isinstance(args.stride, int):
        assert 6 <= args.stride <= 48, "Invalid stride = {} -> 6 <= stride <= 48".format(args.stride)
        data_params["stride"] = args.stride

    if args.gpuid is not None:
        test_params["gpu_id"] = args.gpuid
    else:
        test_params["gpu_id"] = "auto"

    if args.batchsize is not None:
        test_params["batch_size"] = args.batchsize

    if args.usecpu:
        test_params["use_gpu"] = False
        test_params["gpu_id"] = ""

    if args.pyinterp:
        data_params["interp"] = "py"

    # Run stage 1
    print("# Start inference stage1", flush=True)
    inference_stage1(dir_map, contour, os.path.join(dir_models, "modelsx"), dir_out, data_params, test_params)
    print("# Done  inference stage1", flush=True)


    # Run stage 2
    print("# Start inference stage2", flush=True)
    inference_stage2(os.path.join(dir_out, "na.mrc"), contour, os.path.join(dir_models, "modelsx"), dir_out, data_params, test_params)
    print("# Done  inference stage2", flush=True)

    # End
    end = time.time()
    print("# Time consuming {:.4f}".format(end - start), flush=True)
