#!/bin/bash
# Copyright (C) 2024 Tao Li, et al. Huazhong University of Science and Technology

# Users need to properly set the following variables after installation
#######################################################################
    # Python environment of EM2NA
    activate=""
    EM2NA_env=""
    LKH_dir=""
    EM2NA_home=""
#######################################################################
if [ "x${activate}" = "x" ]; then
    echo "Please set the activate first!"
    exit 1
fi
if [ "x${EM2NA_env}" = "x" ]; then
    echo "Please set the EM2NA_env first!"
    exit 1
fi
if [ "x${LKH_dir}" = "x" ]; then
    echo "Please set the LKH_dir first!"
    exit 1
fi
if [ "x${EM2NA_home}" = "x" ];then
    echo "Please set the EM2NA_home first!"
    exit 1
fi

# Get current dir
base_dir=${EM2NA_home}
echo "EM2NA base directory: ${base_dir}"
lkh_dir=${LKH_dir}
echo "LKH locate at ${lkh_dir}"
echo ""

version="1.3"
# Help
print_help(){
    echo "EM2NA v${version} - Automated DNA/RNA modeling from cryo-EM maps"
    echo "Tao Li, Sheng-You Huang et al. 2024"
    echo ""
    echo "Usage: /path/to/EM2NA_v${version}/EM2NA.sh input_map.mrc output_dir [Options]"
    echo "Description:"
    echo "       input_map.mrc : Input cryo-EM density map in MRC2014 format"
    echo "       output_dir    : Output directory"
    echo "                       WARNING if you launch >1 job(s) at a time"
    echo "                       The output_dir **MUST** be different for each job"
    echo "Options:"
    echo "       --seq : Path to input sequence(s) in .fasta format"
    echo ""
    echo "       --contour : Contour level of input map, voxels below will be ignored"
    echo "                   default: '1e-6'"
    echo ""
    echo "       -g : ID(s) of GPU devices to use.  e.g. '0' for GPU #0, and '2,3,6' for GPUs #2, #3, and #6"
    echo "            default: '0'"
    echo ""
    echo "       -b : Number of boxes input into EM2NA in one batch. Users can adjust 'batch_size' according to the VRAM of their GPU devices"
    echo "            Empirically, a GPU with 40 GB VRAM can afford a 'batch_size' of 160"
    echo "            default: '40'"
    echo ""
    echo "       -s : Map grid chunking stride"
    echo "            default: '12'"
    echo ""
    echo "       --usecpu : Specify to run EM2NA deep learning prediction on CPU instead of GPU"
    echo "                  Strongly not recommended since it's very slow"
    echo ""
    echo "       --ncpu : Number of cpus to use to accelarate traces threading"
    echo "                default: '4'"
    echo ""
    echo "       --natype : Nucleic-acid type ['DNA', 'RNA', or 'AUTO'], if 'AUTO', automatically detected by program"
    echo "                  default: 'AUTO'"
    echo ""
    echo "       --time_limit : Running time limit (in seconds) for atomic modeling stage"
    echo "                      default: '172800'"
    echo ""
    echo "       --pyinterp : Interpolation using Python not Fortran"
    echo ""
    #echo "       --help : Specify to see all options"
    #echo ""
    echo "       --keep_temp_files : Keep intermediate files"
    echo "                           Specify to keep temp files"
    echo ""
    #echo ""
    #echo "       --skip_dl : Skip deep learning"
    #echo "                   Specify to skip deep learning detection"
    #echo ""
    #echo "       --skip_model : Skip modeling"
    #echo "                      Specify to only do deep learning detection"
    #echo ""
}

# Print help
if [ $# -lt 2 ];then
    print_help
    exit 1
fi


#############################
# Start program #############
#############################
echo "Job start at" `date`


# Activate env if needed
. ${activate} ${EM2NA_env}

# Read input
input_map=$1
output_dir=$2

# Parameters
input_seq=None
input_ss=None
contour=1e-6

bsize=40
gpuid="0"
ncpu=4
usecpu=""
natype="AUTO"
mode="multimer"
skip_dl="False"
skip_model="False"
keep_temp_files="True"
stride=12
time_limit="172800"
pyinterp=""
use_pyinterp="False"

# Read user-parameters
while [ $# -ge 3 ];do
    case $3 in
    --seq)
        shift
        input_seq=$3;;
    --ss)
        shift
        input_ss=$3;;
    --contour)
        shift
        contour=$3;;
    -g)
        shift
        gpuid=$3;;
    -b)
        shift
        bsize=$3;;
    -s)
        shift
        stride=$3;;
    --natype)
        shift
        natype=$3;;
    --ncpu)
        shift
        ncpu=$3;;
    --usecpu)
        usecpu="--usecpu";;
    --keep_temp_files)
        keep_temp_files="True";;
    --skip_dl)
        skip_dl="True";;
    --skip_model)
        skip_model="True";;
    --pyinterp)
        pyinterp="--pyinterp"
        use_pyinterp="True";;
    --time_limit)
        shift
        time_limit=$3;;
    *)
        echo " ERROR: wrong command argument \"$3\" !!"
        echo " Type \"$0\" for help !!"
        exit 2;;
    esac
    shift
done

temp_dir=${output_dir}/temp
echo "All intermediate results (including logs) will be written to: ${temp_dir}"
if [ ! -e ${output_dir} ]; then
    mkdir ${output_dir}
fi
if [ ! -e ${temp_dir} ]; then
    mkdir ${temp_dir}
fi


# Check input
echo "Please check the input"
echo "Input map -> ${input_map}"
echo "Input map contour -> ${contour}"
echo "Input sequence  -> ${input_seq}"
echo "Map contour -> ${contour}"
echo "Run on GPU -> ${gpuid}"
echo "Batch size -> ${bsize}"
echo "Stride -> ${stride}"
echo "NA type -> ${natype}"
echo "Skip dl -> ${skip_dl}"
echo "Skip model -> ${skip_model}"
echo "Keep temp files -> ${keep_temp_files}"
echo "Pyinterp -> ${use_pyinterp}"
echo "Time limit ->" ${time_limit}
sleep 1s
echo "" > ${temp_dir}/temp.log


# Deep learning detection
if [ ${skip_dl} = "False" ]; then
    echo "Deep-learning prediction will take several minutes for a map of ~(300 Angstrom)^3"
    # Deep learning detection
    echo "Start stage 1 - deep learning detection"
    python ${base_dir}/preds.py -c ${contour} -i ${input_map} -o ${temp_dir} -s ${stride} -b ${bsize} -m ${base_dir} -g ${gpuid} ${usecpu} ${pyinterp}
    if [ -e ${temp_dir}/p.mrc ]; then
        echo "Deep learning detection complete!"
    else
        echo "Deep learning detection failed, please check inputs"
        echo "If Error occurs with correct inputs, please report bugs!"
        exit 1
    fi
fi


# Modeling
if [ ${skip_model} = "False" ]; then
    echo "Start modeling"
    echo "Structure modeling will take < 5 minute for < 500  nucleotides"
    echo "                             ~ 1 hour   for ~ 1500 nucleotides"
    echo "                             + more time for larger structures"

    # Convert to points
    ${base_dir}/bin/getp --in ${temp_dir}/p.mrc  --out ${temp_dir}/p.pdb  --thresh 3 --nt ${ncpu} >> ${temp_dir}/temp.log
    ${base_dir}/bin/getp --in ${temp_dir}/c4.mrc --out ${temp_dir}/c4.pdb --thresh 3 --nt ${ncpu} >> ${temp_dir}/temp.log
    ${base_dir}/bin/getp --in ${temp_dir}/n.mrc  --out ${temp_dir}/n.pdb  --thresh 3 --nt ${ncpu} >> ${temp_dir}/temp.log

    if [ -s ${temp_dir}/p.pdb ]; then
        echo "Find atoms complete!"
    else
        echo "Find atoms failed, please report bugs!"
        exit 1
    fi

    python ${base_dir}/build.py -lkh ${lkh_dir} -lib ${base_dir} -seq ${input_seq} -p ${temp_dir}/p.pdb -c ${temp_dir}/c4.pdb -n ${temp_dir}/n.pdb -pmap ${temp_dir}/p.mrc -cmap ${temp_dir}/c4.mrc -nmap ${temp_dir}/n.mrc -amap ${temp_dir}/aa.mrc -o ${temp_dir} -natype ${natype} ${pyinterp} --time_limit ${time_limit} >> ${temp_dir}/temp.log

    if [ ${input_seq} != "None" ]; then
        echo "Refine seq assign"
        # Refine seq. assign
        python ${base_dir}/opt_seq_assign.py -p ${temp_dir}/denovo.cif --aamap ${temp_dir}/aa.mrc --aaprob ${temp_dir}/prob.npz --seq ${input_seq} -l ${base_dir} -o ${temp_dir} ${pyinterp} >> ${temp_dir}/temp.log

        # Finalize
        cp ${temp_dir}/reassign_helix.cif ${temp_dir}/output.cif

    else
        cp ${temp_dir}/denovo_without_seq.cif ${temp_dir}/output.cif
    fi

    # Score
    python ${base_dir}/eval_local.py -p ${temp_dir}/output.cif -pmap ${temp_dir}/p.mrc -cmap ${temp_dir}/c4.mrc -nmap ${temp_dir}/n.mrc -o ${temp_dir}/scored.cif ${pyinterp}

    # Final
    cp ${temp_dir}/scored.cif ${output_dir}/output.cif

    # Remove temp files
    if [ ${keep_temp_files} = "False" ]; then
        echo "Clean temp files"
        rm -rf ${temp_dir}
    else
        echo "Keep temporary files at ${temp_dir}"
    fi

    echo "Complete, please check output model named output.cif at ${output_dir}"
fi


########################
# End program #########
########################
echo "Job end   at" `date`
