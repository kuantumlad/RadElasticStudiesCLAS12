import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import exclusive.ENPip

def hmm2 = new H1F("hmm2", "missing mass squared", 200,0,10)

def beam = new Particle(11, 0,0,10.6)//7.546)
def target = new Particle(2212, 0,0,0)

for(fname in args) {
def reader = new HipoDataSource()
reader.open(fname)

while(reader.hasEvent()) {
  def event = reader.getNextEvent()
  if (event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter")) {
    def (ele, neu, pip) = ENPip.getENPip(event)*.particle
    if(ele!=null) {
      def eX = new Particle(beam)
      eX.combine(target, 1)
      eX.combine(ele,-1)
      eX.combine(pip,-1)

      hmm2.fill(eX.mass2())
    }
  }
}

reader.close()
}

def out = new TDirectory()
out.mkdir('/enpip')
out.cd('/enpip')
out.addDataSet(hmm2)
out.writeFile('enpip_out.hipo')
