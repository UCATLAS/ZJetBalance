# setupATLAS
# localSetupPyAMI


import pyAMI.client
import pyAMI.atlas.api as AtlasAPI


class MyDxAODInfoGetterOk():
    
    ################################################
    def __init__(s):
        s.client = pyAMI.client.Client('atlas')
        AtlasAPI.init()
    
    ################################################
    def getInfo(s,dsid,daodType,ptag,rtag=None):
        daod=s.getLdn(dsid,daodType,ptag,rtag)
        
        if not daod:
            print 'ERROR> in getLdn, please find your query is enough to select one DS'
            print '       Currently follows are used: dsid=',dsid,'daodType=',daodType,'ptag=',ptag,'rtag=',rtag
            return
        
        aod  = s.getOriginalAOD(daod)
        
        if not aod:
            print 'ERROR> no AOD File is not found for ldn=',ldn
            return
        
        aodInfo  = s.encodeDSInfo(aod)
        daodInfo = s.encodeDSInfo(daod)
        
        aodEvents  = float(aodInfo['events'])
        daodEvents = float(daodInfo['events'])
        skimmingRate = daodEvents/aodEvents
        return {'events':daodInfo['events'], 'xsec':daodInfo['xsec'], 'fit_eff':daodInfo['fit_eff'], 'daod_ldn':daod, 'aod_ldn':aod, 'dsid':dsid, 'skim_eff':skimmingRate, 'daod_events':daodEvents, 'aod_events':aodEvents}
    
    ################################################
    def encodeDSInfo(s, ldn):
        datasetinfo      = AtlasAPI.get_dataset_info(s.client,ldn)
        neventsAMI       = float(datasetinfo[0]['totalEvents'])
        crosssectionAMI  = float(datasetinfo[0]['crossSection'])
        if 'approx_GenFiltEff' in datasetinfo[0].keys():
            filtereffAMI     = float(datasetinfo[0]['approx_GenFiltEff'])
        else :
            filtereffAMI     = None
        return {'events':neventsAMI, 'xsec':crosssectionAMI, 'fit_eff':filtereffAMI, 'ldn':ldn}
    
    ################################################
    def getLdn(s,dsid,daodType,ptag,rtag=None):
        pattern=''
        if not rtag:
            pattern='mc15_13TeV.'+dsid+'%.merge.DAOD_'+daodType+'%'+ptag+'%'
        else :
            pattern='mc15_13TeV.'+dsid+'%.merge.DAOD_'+daodType+'%'+rtag+'%'+ptag+'%'
            
        samples=AtlasAPI.list_datasets(s.client, patterns = [pattern])
        
        if len(samples)==0:
            print 'NO DAOD_%s DS found for %s'%(DAOD_Type, dsid)            
            return None
        elif len(samples)!=1:
            print 'More than one DAOD_%s DS round for %s'%(DAOD_Type, dsid)
            for sample in samples:
                print sample
            print 'please set rtag to select one proper one'
            return None
        
        return samples[0]['ldn']
    
    ################################################
    def getOriginalEVNT(s,ldn):
        provs=AtlasAPI.get_dataset_prov(s.client, ldn)
        
        evnt=None
        
        for sample in provs['node']:
            if 'evgen.EVNT' in sample['logicalDatasetName']:
                evnt=sample['logicalDatasetName']
        
        return evnt
    
    ################################################
    def getOriginalAOD(s,ldn):
        provs=AtlasAPI.get_dataset_prov(s.client, ldn)
        
        evnt=None
        
        for sample in provs['node']:
            if '.AOD.' in sample['logicalDatasetName']:
                evnt=sample['logicalDatasetName']
                break;
            
        return evnt
    
if __name__ == '__main__':
    getter = MyDxAODInfoGetterOk()
    info=getter.getInfo(dsid='361372',daodType='SUSY2',ptag='p2375')
