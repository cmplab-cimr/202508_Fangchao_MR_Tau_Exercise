# 202508_Fangchao_MR_Tau_Exercise
Flux analysis of fruit flies for Fangchao to study effects of methionine restriction, taurine and exercise.


## Getting started

This script could also be executed as a raw Python project. Ensure that you have the latest version of Anaconda installed (Python 3.12 or higher). The script can be run in the `base` environment, but you will also need to install the `xlsxwriter` package. First switch to a target directory and download the source code:

```shell script
git clone https://github.com/cmplab-cimr/202508_Fangchao_MR_Tau_Exercise
```

Switch to the source direct, add PYTHONPATH environment and run the `main.py`:

```shell script
cd 202508_Fangchao_MR_Tau_Exercise
export PYTHONPATH=$PYTHONPATH:`pwd`
python main.py
```

You could try multiple different arguments according to help information. For example:


```shell script
python main.py computation experimental_mfa flux_analysis -t
```
