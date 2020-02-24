package exclusive

import org.jlab.io.hipo.HipoDataEvent
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import pid.electron.Electron
import pid.proton.Proton

class EP {
  static def getEP(HipoDataEvent event) {
    def partbank = event.getBank("REC::Particle")
    def calbank = event.getBank("REC::Calorimeter")

    def eleind = Electron.findElectron(event)
    if(eleind!=null) {
      //find proton with maximum energy
      def proind = Proton.findProton(event)

      def secs = [calbank.getShort('pindex')*.toInteger(), calbank.getByte('sector')].transpose().collectEntries()

      if(proind!=null) {
        def electron = new Particle(11, *['px', 'py', 'pz'].collect{partbank.getFloat(it, eleind)})
        def proton = new Particle(2212, *['px', 'py', 'pz'].collect{partbank.getFloat(it, proind)})

        def status = partbank.getShort('status', proind)
        return [[particle:electron, pindex:eleind, sector:secs[eleind]],
          [particle:proton, pindex:proind, status:status]]
      }
    }
    return [null, null]
  }


  static def getElastic(HipoDataEvent event) {
    def (electron, proton) = getEP(event)
    if(electron!=null) {
      //angle between electron and proton planes (should be close to 180)
      def zaxis = new Vector3(0,0,1)
      def enorm = electron.particle.vector().vect().cross(zaxis)
      def pnorm = proton.particle.vector().vect().cross(zaxis)
      def phi = enorm.theta(pnorm)
      if(phi>176) return [electron, proton]
    }
    return [null, null]
  }
}

