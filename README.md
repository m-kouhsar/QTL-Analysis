# eQTL Analysis using matrixEQTL R package
Run the first step using one of the following commands:
```
# Submit as a slurm job:
sbatch ./Slurm/1_MatrixEQTL.PrepareData_Slurm.sh ./MatrixEQTL.config

# Run as a simple bash:
bash ./1_MatrixEQTL.PrepareData.sh ./MatrixEQTL.config
```
You can run second and third stepd in the same way. 

There are some example input data in the [Example](https://github.com/m-kouhsar/QTL-Analysis/tree/main/Example) folder. Use the following commands to run the pipeline on the example data:
```
git clone https://github.com/m-kouhsar/QTL-Analysis.git
cd QTL-Analysis/
bash 1_MatrixEQTL.PrepareData.sh MatrixEQTL.config
bash 2_MatrixEQTL.Main.sh MatrixEQTL.config
bash 3_MatrixEQTL.SaveResults.sh MatrixEQTL.config
```
