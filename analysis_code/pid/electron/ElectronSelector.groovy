package pid.electron
import event.Event
import pid.electron.ElectronFromEvent

class ElectronSelector {

    def electronCuts = new ElectronFromEvent()
    def electronCutStrategies

    ElectronSelector(){
	this.initializeCuts()
    }
    
    def initializeCuts(){
	this.electronCutStrategies = [
	    electronCuts.passElectronStatus,
	    electronCuts.passElectronChargeCut,
	    electronCuts.passElectronEBPIDCut,
	    electronCuts.passElectronNpheCut,
	    electronCuts.passElectronVertexCut,
	    electronCuts.passElectronPCALFiducialCut,
	    electronCuts.passElectronEIEOCut,
	    electronCuts.passElectronDCR1,
	    electronCuts.passElectronDCR2,
	    electronCuts.passElectronDCR3
	]
    }

    def getGoodElectronWithCuts(event){
	//return map - key is index of track in REC::Particle and value is list of booleans for the cuts
	def el_cut_result = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, electronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()	  		
	return el_cut_result
    }

    def getGoodElectron(event){
	//return a list of REC::Particle indices for tracks passing all electron cuts
	def el_cut_result = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, electronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()
	return el_cut_result.findResults{el_indx, cut_result -> !cut_result.contains(false) ? el_indx : null}
    }

}
