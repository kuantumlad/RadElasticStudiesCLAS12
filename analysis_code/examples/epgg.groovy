#!/home/kenjo/.groovy/coatjava/bin/run-groovy

import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import exclusive.EPGG

def hpi0m = new H1F("hpi0m", "pi0 mass", 200,0,0.3)
def hmm2 = new H1F("hmm2", "missing mass squared", 200,-2,4)

def beam = new Particle(11, 0,0,10.6)//7.546)
def target = new Particle(2212, 0,0,0)

for(fname in args) {
def reader = new HipoDataSource()
reader.open(fname)

while(reader.hasEvent()) {
  def event = reader.getNextEvent()
  if (event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter")) {
    def (ele, pro, g1, g2) = EPGG.getEPGG(event)*.particle
    if(ele!=null) {
      def pi0 = new Particle(g1)
      pi0.combine(g2,1)
      hpi0m.fill(pi0.mass())

      def epX = new Particle(beam)
      epX.combine(target, 1)
      epX.combine(ele, -1)
      epX.combine(pro, -1)
      hmm2.fill(epX.mass2())
    }
  }
}

reader.close()
}

def out = new TDirectory()
out.mkdir('/epgg')
out.cd('/epgg')
out.addDataSet(hpi0m)
out.addDataSet(hmm2)
out.writeFile('epgg_out.hipo')
