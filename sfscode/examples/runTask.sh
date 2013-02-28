#!sh
#$ -S /bin/bash
#$ -t 1-10
#$ -cwd
TASKFILE=tasklist
TASK=$(cat $TASKFILE | head -n $SGE_TASK_ID | tail -n 1)
echo $TASK
$TASK -o out_$SGE_TASK_ID.txt

exit 0
