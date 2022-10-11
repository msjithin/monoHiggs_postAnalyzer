set -e

cd postAnalyzer

echo "hadding files......"
sh hadd_files.sh

cd plotting_script

echo "Lumi-scaling and arranging histograms...."
bash execute_main.sh


