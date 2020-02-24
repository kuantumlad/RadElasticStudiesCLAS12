#!/home/kenjo/.groovy/coatjava/bin/run-groovy

import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import exclusive.EPG

def hmm2 = new H1F("hmm2", "missing mass squared", 200,-2,4)

def beam = new Particle(11, 0,0,10.6)//7.546)
def target = new Particle(2212, 0,0,0)

for(fname in args) {
def reader = new HipoDataSource()
reader.open(fname)

while(reader.hasEvent()) {
  def event = reader.getNextEvent()
  if (event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter")) {
    def (ele, pro, gam) = EPG.getEPG(event)*.particle
    if(ele!=null) {
      def epX = new Particle(beam)
      epX.combine(target, 1)
      epX.combine(ele,-1)
      epX.combine(pro,-1)

      hmm2.fill(epX.mass2())
    }
  }
}

reader.close()
}

def out = new TDirectory()
out.mkdir('/epg')
out.cd('/epg')
out.addDataSet(hmm2)
out.writeFile('epg_out.hipo')
