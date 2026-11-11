## DNN with 3 features - Train, Evaluation, Draw Histograms
### 0. Presets
```{.Bash}
  source /cvmfs/sft.cern.ch/lcg/views/LCG_108_cuda/x86_64-el9-gcc13-opt/setup.sh
  module load cuda/12.5
  module load cudnn/9.3
```
 
### 1. Train
#### 1.1. Train 2 separate binary class classification 
```{.Bash}
python train.py
```
**Train** is performed using processed ntuples in `../nanoaodframe_*LFV/` and corresponding part should be **MODIFIED**.

Training is performed with signal (LFV sample mixture) and background (tt semileptonic and dileptonic mixture) samples.

When training is done, `./<train output folder>/` will be generated (e.g. TTnov_01_norm_20nodes_3layers_s1b1) with result plots (correlation, loss, accuracy, DNN output, confusion matrix).

Best model is stored in the output folder as `best_model.h5`.
#### 1.2. Train 2 separate multi class classification 
```{.Bash}
python train_multi.py
```


### 2. Evaluation, post processing, drawing
#### 2.1. Multi class
```{.Bash}
python eval_multi.py -I top_lfv_multiClass_March2024_AfterPreAppTalk -O DNN_out_0424
```
On top of evaluated folder, we run usual post processing scripts and them collect all histograms for combine tool
```{.Bash}
python ../nanoaodframe/postprocess.py -I DNN_out_0424 -Y 2018
#If you produce fine binned histograms for smoothing but want to draw histos in the original bin:
python ../nanoaodframe/postprocess.py -I DNN_out_0424 -Y 2018 --postfix _orgbin

python ../nanoaodframe/print_syst_table.py -I DNN_out_0424 --postfix _orgbin -D -Y 2018

cd DNN_out_0424/2018_postprocess_orgbin
../../../plotIt/plotIt -o ../figure_2018/ ../../../plotIt/configs/TOP-22-011/config_2018.yml -y -s
../../../plotIt/plotIt -o ../figure_2017/ ../../../plotIt/configs/TOP-22-011/config_2017.yml -y -s
../../../plotIt/plotIt -o ../figure_2016post/ ../../../plotIt/configs/TOP-22-011/config_2016post.yml -y -s
../../../plotIt/plotIt -o ../figure_2016pre/ ../../../plotIt/configs/TOP-22-011/config_2016pre.yml -y -s

python postprocess_2.py -I DNN_out_0424
python ../nanoaodframe/stack_signals_v2.py -I DNN_out_0424 --postfix _orgbin
python ../nanoaodframe/plot_run2.py -I DNN_out_0424 --postfix _orgbin -D
```
**Evaluation** is performed with the best model from training and produces DNN output score histograms (`h_dnn_pred_S5` and all systematic variations) with data and MC samples.

To use the proper model, check the directory of the processed ntuple used in training. The 24 Run 3 input features per channel (`muon` / `electron`) are defined in `utils/feature_config.py`.

The output folder is defined in `eval_multi.py` (`eval_multi_bi.py`). During evaluation, `hcounter`, `ScaleWeightSum`, `PSWeightSum`, and `LHEPdfWeightSum` are automatically copied from the input process ROOT files so that `postprocess.py` can renormalize the theory scale (`mescale`, `renscale`, `facscale`), PS (`isr`, `fsr`), and PDF/$\alpha_S$ variations by their pre-selection generator weight sums and merge external shape variations into a single nominal ROOT file.

The plotIt does not support stack for signals, thus it is done by dedicated code.

The `plot_run2.py` will read information from existing configs for each year, then combine into one, usint templates.

### Advanced
> To modify some details of the plots from training, cross section for normalization, histogram styles, modules in `utils/` folder will be helpful.
