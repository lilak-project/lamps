f1 = open("../common/LAMPS_TPC_PAD_info_11");
for ulayer in range(0,42):#[37,38,39,40,41]:
    f1.seek(0)
    countChannel = 0
    while True:
        line = f1.readline()
        line = line.strip()
        if not line or len(line)==0:
            break
        sline = line.split()
        #print(sline)
        ch = int(sline[3])
        x = float(sline[4])
        y = float(sline[5])
        section = int(sline[6])
        padid = int(sline[7])
        layer = int(sline[8])
        if x>0 and section==1 and layer==ulayer:
            countChannel = countChannel + 1
    print(ulayer, countChannel)
