import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import pid.electron.ElectronFromEvent
import event.Event
import event.EventConverter
import utils.KinTool
import pid.electron.ElectronSelector

//use these two lines for the second and third method
def ele_selector = new ElectronSelector()
//ele_selector.initializeCuts()

//use two lines below for first method
def electron = new ElectronFromEvent();
def myElectronCutStrategies = [
    electron.passElectronStatus,
    electron.passElectronChargeCut,
    electron.passElectronEBPIDCut,
    electron.passElectronNpheCut,
    electron.passElectronVertexCut,
    electron.passElectronPCALFiducialCut,
    electron.passElectronEIEOCut,
    electron.passElectronDCR1,
    electron.passElectronDCR2,
    electron.passElectronDCR3
]
   
for(fname in args) {
    def reader = new HipoDataSource()
    reader.open(fname)
    
    while(reader.hasEvent()) {
	def data_event = reader.getNextEvent()
	def event = EventConverter.convert(data_event)

	// first method for selecting electrons illustrates what the ElectronSelector class is doing under the hood.
	def my_el_cuts = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()	  
	my_el_cuts.each{ index, value ->
	    def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
	    def p = lv.mag()
	    def vz = event.vz[index]
	    def theta = Math.toDegrees(lv.theta())
	    def phi = Math.toDegrees(lv.phi())
	    //do other stuff here
	}

	// second method here will return a list of indicies for tracks passing all electron cuts.
	// decide which of the electron candidates you want here
	def my_good_el = ele_selector.getGoodElectron(event)	
	// third method will return a map with key as the REC::Particle index and value as the list of booleans describing if the track passed the cut or not.
	def my_good_el_with_cuts = ele_selector.getGoodElectronWithCuts(event)

    }
}

reader.close()
