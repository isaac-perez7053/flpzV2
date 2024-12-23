set modeName = $1


set numJobs_ind = `grep NumSAM headerFile_$modeName.dat | awk '{print $2}'`
#set numJobs = 5 #how many jobs to split into, remember indexing


mkdir SMODES_$modeName
mv    headerFile_$modeName.dat SMODES_$modeName
#set numJobs_ind=`echo $numJobs\- 1 | bc`
# rm joblist
foreach II (`seq 0 $numJobs_ind` )
        cp -r boilerplate SMODES_$modeName/dist_${II}
        sed -i "/CELLDEF/ r dist_${modeName}_${II}" SMODES_${modeName}/dist_${II}/template.abi 
        sed -i "s/CELLDEF/ /g" SMODES_$modeName/dist_${II}/template.abi 
	rm dist_${modeName}_${II}
	mv SMODES_$modeName/dist_${II}/template.abi SMODES_$modeName/dist_${II}/dist_${II}.abi
	sed -i "s/DISTNAME/dist_${II}/g" SMODES_$modeName/dist_${II}/jobscript.sh
        echo "cd SMODES_$modeName/dist_${II}; sbatch jobscript.sh; cd -" >> joblist
end

