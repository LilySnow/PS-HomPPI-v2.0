#for i in `cat /home/lixue/PSHOMPPIv1.3/uploadData/BM5_dimers/BM5_dimer.lst `;do
#for i in `cat /home/lixue/PSHOMPPIv1.3/uploadData/BM5_dimers/huge_case1 `;do
for i in `egrep -v '^#' /data/web_servers/lxue/PSHOMPPIv2.0/uploadData/BM5_dimers/BM5_dimer_hugeCasesExcl.lst `;do
#qsub -q short -V <<END
#    cd /home/lixue/PSHOMPPIv1.3/WEB-INF/cgi_BM5_Dec22nd2016
    perl PSHomPPI.cgi $i

#END
sleep 30
done


