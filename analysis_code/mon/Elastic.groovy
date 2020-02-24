package mon
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H2F
import exclusive.EP

class Elastic {
  def beam = new Particle(11, 0,0,10.6)//7.546)
  def target = new Particle(2212, 0,0,0)

  def plots = [:].withDefault{[:].withDefault{new H2F("hemm2fi_$it", "MM2 vs phi", 250,140,190, 250,0,5)}}
  def data = [:].withDefault{[:].withDefault{0}}

  def dets = ['FT', 'FD', 'CD']

  def processEvent(event) {
    if (event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter")) {
      def (ele, pro) = EP.getEP(event)
      if(ele!=null) {
        def eX = new Particle(beam)
        eX.combine(target, 1)
        eX.combine(ele.particle, -1)

        //angle between electron and proton planes (should be close to 180)
        def zaxis = new Vector3(0,0,1)
        def enorm = ele.particle.vector().vect().cross(zaxis)
        def pnorm = pro.particle.vector().vect().cross(zaxis)
        def phi = enorm.theta(pnorm)

        def calbank = event.getBank("REC::Calorimeter")
        def elesec = (0..<calbank.rows()).collect{
          (calbank.getShort('pindex',it).toInteger() == ele.pindex &&
          calbank.getByte('detector',it).toInteger() == DetectorType.ECAL.getDetectorId()) ? calbank.getByte('sector',it) : null
        }.find()
        int prost = pro.status/2000
        /**/

        if(elesec!=null) {
          def id = 'pro@'+dets[prost]
          plots[id]["before_sec$elesec"].fill(phi, eX.mass2())
          if(phi>176 && eX.mass2()<1.5) {
            data[id]["sec$elesec"]++
            plots[id]["after_sec$elesec"].fill(phi, eX.mass2())
          }
        }
      }
    }
  }
}
