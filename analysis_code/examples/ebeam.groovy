import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.TDirectory
import exclusive.EP

def hbeam = [:].withDefault{new H1F("hbeam_$it", "beam energy for $it", 200,5,15)}
def heth = new H1F("heth", "electron theta", 250, 0, 25)

def beam = new Particle(11, 0,0,10.6)//7.546)
def target = new Particle(2212, 0,0,0)

for(fname in args) {
def reader = new HipoDataSource()
reader.open(fname)

while(reader.hasEvent()) {
  def event = reader.getNextEvent()
  if (event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter")) {
    def (ele, pro) = EP.getEP(event)
    if(ele!=null) {
      //angle between electron and proton planes (should be close to 180)
      def zaxis = new Vector3(0,0,1)
      def enorm = ele.particle.vector().vect().cross(zaxis)
      def pnorm = pro.particle.vector().vect().cross(zaxis)
      def phi = enorm.theta(pnorm)

      def ebeam = 0.938*(1/(Math.tan(ele.particle.theta()/2)*Math.tan(pro.particle.theta()))-1)
      if(phi>176) {
        def theta = Math.toDegrees(ele.particle.theta())
        heth.fill(theta)
        int imom = (heth.getXaxis().getBin(theta)-1)/50
        hbeam["sec${ele.sector}_mom${imom}"].fill(ebeam)
      }
    }
  }
}

reader.close()

def out = new TDirectory()
out.mkdir('/ep')
out.cd('/ep')
hbeam.values().each{out.addDataSet(it)}
out.addDataSet(heth)
out.writeFile('ebeam_out.hipo')

}
