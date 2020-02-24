package pid.proton

import org.jlab.clas.physics.Vector3

class Proton {
  static def findProton = { event ->
    def pbank = event.getBank("REC::Particle")
    return (0..<pbank.rows()).findAll{pbank.getInt('pid',it)==2212}
      .max{ind -> (new Vector3(*['px', 'py', 'pz'].collect{pbank.getFloat(it,ind)})).mag2()}
  }
}
