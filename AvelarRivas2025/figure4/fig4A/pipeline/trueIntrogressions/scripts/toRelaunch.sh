stDir="st/"
mainDir=$(pwd)

mkdir -p data/ascTrees

echo "echo Launching New Trees!" > scripts/launchNewtrees.sh
echo "module load iq-tree/2.3.6" >> scripts/launchNewtrees.sh
for errfile in $(grep "Invalid" ${stDir}*.ne | cut -d":" -f1); do
  echo $errfile
  commandst=$(grep "iqtree2" $errfile | head -1)
  phyFile=$(grep "input" $errfile | cut -d"/" -f4| head -1)
  strain=$(grep "sample" $errfile | cut -d"=" -f2 | cut -d"," -f1 | head -1)
  mkdir -p data/ascTrees/$strain
  cp $mainDir/data/phyl/$strain/$phyFile.varsites.phy data/ascTrees/$strain
  echo $commandst > tmpf
  sed -i "s/phyl/ascTrees/g" tmpf
  sed -i "s/.phy/.phy.varsites.phy/g" tmpf
  cattmpf=$(cat tmpf)
  toMakeTree=$(echo $cattmpf" --redo")
  echo $toMakeTree" > "$stDir$strain"_"${phyFile}".asc.out 2> "$stDir$strain"_"${phyFile}".asc.err" >> scripts/launchNewtrees.sh
  echo "echo 'Done with '${phyFile}" >> scripts/launchNewtrees.sh
done

bash scripts/launchNewtrees.sh
