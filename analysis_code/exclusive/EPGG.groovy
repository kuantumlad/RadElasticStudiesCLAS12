package exclusive

import org.jlab.io.hipo.HipoDataEvent
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3

class EPGG {
  static def getEPPi0(HipoDataEvent event) {
    def partbank = event.getBank("REC::Particle")
    def calbank = event.getBank("REC::Calorimeter")

    def findElectron = { pbank -> (0..<pbank.rows()).find{pbank.getInt('pid',it)==11 && pbank.getShort('status',it)<0} }
    def findProton = { pbank -> (0..<pbank.rows()).findAll{pbank.getInt('pid',it)==2212}
      .max{ind -> (new Vector3(*['px', 'py', 'pz'].collect{pbank.getFloat(it,ind)})).mag2()}
    }
    def findGammas = { pbank -> (0..<pbank.rows()).findAll{pbank.getInt('pid',it)==22}
      .sort{ind -> (new Vector3(*['px', 'py', 'pz'].collect{pbank.getFloat(it,ind)})).mag2()}
      .reverse()
    }

    def inds = []
    for(def findPart in [findElectron, findProton]) {
      def ind = findPart(partbank)
      inds.add(ind)
      if(ind==null)
        return [null, null, null, null]
    }

    def vele = new Particle(11, *['px', 'py', 'pz'].collect{partbank.getFloat(it,inds[0])})
    def gammas = findGammas(partbank).findAll{ind -> vele.vector().vect().theta(new Vector3(*['px', 'py', 'pz'].collect{partbank.getFloat(it,ind)}))>5}
    def pi0s = gammas.dropRight(1).withIndex().collect{ig1, ind1->
      gammas.drop(ind1+1).findResults{ig2->
        def vg1 = new Particle(22, *['px','py','pz'].collect{partbank.getFloat(it,ig1)})
        def vg2 = new Particle(22, *['px','py','pz'].collect{partbank.getFloat(it,ig2)})
        def vpi0 = new Particle(vg1)
        vpi0.combine(vg2,1)

        def isGoodPi0 = vg1.vector().vect().theta(vg2.vector().vect()) > 2
        isGoodPi0 ? [ig1, ig2, vpi0] : null
      }
    }.collectMany{it}.sort{it[2].e()}.reverse()

    if(pi0s.size()<1)
      return [null, null, null, null]
    inds.addAll(pi0s[0][0..1])

    def secs = [calbank.getShort('pindex')*.toInteger(), calbank.getByte('sector')].transpose().collectEntries()

    def parts = [11,2212,22,22].withIndex()
      .collect{pid,i -> new Particle(pid, *['px', 'py', 'pz'].collect{partbank.getFloat(it, inds[i])}) }

    return (0..<4).collect{[particle:parts[it], pindex:inds[it], sector:secs[inds[it]]]}
  }
}

