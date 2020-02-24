package exclusive

import org.jlab.io.hipo.HipoDataEvent
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3

class ENPip {
  static def getENPip(HipoDataEvent event) {
    def partbank = event.getBank("REC::Particle")
    def calbank = event.getBank("REC::Calorimeter")

    def findElectron = { pbank -> (0..<pbank.rows()).find{pbank.getInt('pid',it)==11 && pbank.getShort('status',it)<0} }
    def findNeutron = { pbank -> (0..<pbank.rows()).findAll{pbank.getInt('pid',it)==2112}
      .max{ind -> (new Vector3(*['px', 'py', 'pz'].collect{pbank.getFloat(it,ind)})).mag2()}
    }
    def findPip = { pbank -> (0..<pbank.rows()).findAll{pbank.getInt('pid',it)==211}
      .max{ind -> (new Vector3(*['px', 'py', 'pz'].collect{pbank.getFloat(it,ind)})).mag2()}
    }

    def inds = []
    for(def findPart in [findElectron, findNeutron, findPip]) {
      def ind = findPart(partbank)
      inds.add(ind)
      if(ind==null)
        return [null, null, null]
    }

    def secs = [calbank.getShort('pindex')*.toInteger(), calbank.getByte('sector')].transpose().collectEntries()

    def parts = [11,2112,211].withIndex()
      .collect{pid,i -> new Particle(pid, *['px', 'py', 'pz'].collect{partbank.getFloat(it, inds[i])}) }

    return (0..<3).collect{[particle:parts[it], pindex:inds[it], sector:secs[inds[it]]]}
  }
}

