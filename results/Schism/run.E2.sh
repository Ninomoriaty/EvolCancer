#! /bin/bash
# run schism in step-through mode

baseDir=$(dirname $BASH_SOURCE)

# assuming the configuration file is in E2.yaml 
# in the same directory as run.E2.sh

currentPath=$(grep 'working_dir' $baseDir/E2.yaml | sed 's/working_dir://g')

# check to see if the path specified as working directory
# exists. If it does, proceed with the analysis. Otherwise,
# create a new configuration file where the path has been 
# substituted with the path to this shell script

echo "The specified working directory:  $currentPath"

if [ -d $currentPath ]; then
    echo "The working directory specified in E2.yaml exists."
    echo "proceeding to analysis"
    # run schism in sequential mode
    conf="$baseDir/E2.yaml"
else
    echo "The working directory specified in E2.yaml does not exists."
    echo "Using the path to this shell script as the working directory."

    echo "# adjusted path to working directory" >> $baseDir/E2.path.yaml
    echo "working_dir: $baseDir" >> $baseDir/E2.path.yaml
    echo "" >> $baseDir/E2.path.yaml
    grep -v "working_dir" $baseDir/E2.yaml >> $baseDir/E2.path.yaml
    # run schism in sequential mode
    conf="$baseDir/E2.path.yaml"
fi


# prepare input of hypothesis test
runSchism prepare_for_hypothesis_test -c $conf

# perform hypothesis test
runSchism hypothesis_test -c $conf

# visualize hypothesis test results
runSchism plot_cpov -c $conf

# perform 10 independent runs of  genetic algorithm
# in parallel by calling each command below on a distinct 
# processor

runSchism run_ga --config $conf --mode parallel --runID 1
runSchism run_ga --config $conf --mode parallel --runID 2
runSchism run_ga --config $conf --mode parallel --runID 3
runSchism run_ga --config $conf --mode parallel --runID 4
runSchism run_ga --config $conf --mode parallel --runID 5
runSchism run_ga --config $conf --mode parallel --runID 6
runSchism run_ga --config $conf --mode parallel --runID 7
runSchism run_ga --config $conf --mode parallel --runID 8
runSchism run_ga --config $conf --mode parallel --runID 9
runSchism run_ga --config $conf --mode parallel --runID 10

# or in serial mode on a single processor
runSchism run_ga --config $conf --mode serial

# gather topologies from independent runs of genetic algorithm
# generate summary log file and plots
runSchism summarize_ga_results --config $conf

# generate the consensus of all maximum fitness trees and visualize
runSchism consensus_tree -c $conf
