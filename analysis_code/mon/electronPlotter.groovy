import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import java.util.concurrent.ConcurrentHashMap
import pid.electron.Electron

def beam = 10.604

//call electron cut constructor
electron = new Electron();

def myElectronCutStrategies = [
    electron.passElectronStatus,
    electron.passElectronChargeCut,
    electron.passElectronEBPIDCut,
    electron.passElectronNpheCut,
    electron.passElectronVertexCut,
    electron.passElectronPCALFiducialCut,
    electron.passElectronEIEOCut,
    electron.passElectronDCR1
]

reqBanks = [ 
    part : "REC::Particle",
    ec   : "REC::Calorimeter",
    cc   : "REC::Cherenkov",    
    traj : "REC::Traj",
    trck : "REC::Track"
]

def out = new TDirectory()

def histos = new ConcurrentHashMap()
histoBuilders = [

    ptheta : {title -> new H2F("$title","$title", 900, 0, beam, 900, 0.0, 65.0 ) },
    phitheta : {title -> new H2F("$title","$title", 900, 180, -180, 900, 0.0, 65.0 ) },
    vzphi : { title -> new H2F("$title","$title", 140, 0, 140, 1000, -100, 100 ) },
    vztheta : { title -> new H2F("$title","$title", 1000, -100, 100, 900, 0.0, 65.0 ) },
    hitxy : { title -> new H2F("$title","$title", 900, -450.0, 450.0, 900, -450.0, 450.0 ) },
    ecsf : { title -> new H2F("$title","$title", 500, 0.0, beam, 0.0, 0.5) },
    nphe : { title -> new H1F("$title","$title", 40, 0.0, 40 ) },
    eieo : { title -> new H2F("$title","$title", 1000, 1.0, 1000, 1.0 ) },
    vz   : { title -> new H1F("$title","$title", 400, -40, 40) },
    w    : { title -> new H1F("$title","$title", 500, 0.5, beam+0.3 ) },
    q2w  : { title -> new H2F("$title","$title", 500, 0.5, beam+0.3, 500, 0, 2 * beam ) },  
    passrate : {title -> new H1F("$title","$title", 8, 0.0, 8.0 ) },
    wtheta : { title -> new H2F("$title","$title", 100, 0.0, 2.5, 100, 0.0, 20.0 ) }

]

for(fname in args) {
    def reader = new HipoDataSource()
    reader.open(fname)
    
    def cc = 0
    
    while(reader.hasEvent() && cc < 20000) {
	
	def event = reader.getNextEvent()

	def banks  = reqBanks.collect{name, bank -> [name, event.getBank(bank)] }.collectEntries()
	def banks_pres = reqBanks.every{name, bank -> event.hasBank(bank) }

	
	if (banks_pres) { 
	    cc+=1
	    
	    def my_el_cuts = (0..<partb.rows()).collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(banks,ii) } ] }.collectEntries()


	    println(my_el_cuts)

	

	}
    }

}
