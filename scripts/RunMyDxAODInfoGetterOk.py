# setupATLAS
# localSetupPyAMI
# voms-proxy-init -voms atlas

from MyDxAODInfoGetterOk import MyDxAODInfoGetterOk


daodType='SUSY2'
ptag='p2375'
rtag=None # add if needed
outputfilename='amiDumpInfo.txt'

DSIDs=[
    # Zee
    '361372', '361373', '361374',
    '361375', '361376', '361377',
    '361378', '361379', '361380',
    '361381', '361382', '361383',
    # Zmumu
    '361396', '361397', '361398',
    '361399', '361400', '361401',
    '361402', '361403', '361404',
    '361405', '361406', '361407',
    # Z tau tau
    '361420', '361421', '361422', 
    '361423', '361424', '361425', 
    '361426', '361427', '361428', 
    '361429', '361430', '361431'
    ]

getter = MyDxAODInfoGetterOk()

content1=[]
content2=[]
content3=[]

for dsid in DSIDs :
    print 'information for '+dsid+' is being dumped..'
    info=getter.getInfo(dsid=dsid,daodType=daodType,ptag=ptag,rtag=rtag)
    # for DS list to be used for submission
    content1.append(info['daod_ldn'])
    
    # for XS file
    # 361106 1.9506E+06 1.0000E+00 19993000
    unit_pb=1e6
    content2.append('#'+info['daod_ldn'])
    content2.append('%s %e %e %d'%(info['dsid'], float(info['xsec'])*unit_pb, float(info['fit_eff']), int(info['aod_events'])) )
    
    # used in the final protting code                
    content3.append('dXAODSkimmingEfficiency[%d]=%f'%(int(info['dsid']), float(info['skim_eff'])) )


outputfile=open(outputfilename, 'w')

for line in content1:
    print >>outputfile,line

print >>outputfile, ''
for line in content2:
    print >>outputfile,line

print >>outputfile, ''
for line in content3:
    print >>outputfile,line
