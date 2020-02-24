import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentLinkedQueue

class epkpkn {
  def hists = new ConcurrentHashMap()
  def beam = LorentzVector.withPID(11,0,0,10.6)
  def target = LorentzVector.withPID(2212,0,0,0)
  def hmissemm2 = {new H2F("$it","$it",100,-2,5,100,0,5)}
  def hmissephim = {new H2F("$it","$it",100,-2,5,100,0,5)}
  def hphim = {new H1F("$it","$it",75,0.8,1.5)}
  def hcombs = {new H1F("$it","$it",9,1,10)}
  def hW = {new H1F("$it","$it",100,1,4)}
  def hmissmm2phim = {new H2F("$it","$it",100,-1,1,100,0,5)}
  def hmissmm2 = {new H1F("$it","$it",100,-1,1)}

  def processEvent(event) {
    if(event.hasBank("REC::Particle")) {
      def partb = event.getBank("REC::Particle")
      def combs = (0..<partb.rows()).findAll{partb.getShort('status',it)<0 && partb.getInt('pid',it)==11}
        .findResults{iele->
          def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
          return ele.e() > 1.5 ? ele : null
        }.collectMany{ele->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==2212}.findResults{ipro->
            def pro = LorentzVector.withPID(2212,*['px','py','pz'].collect{partb.getFloat(it,ipro)})
            return (pro.p()>0.5 && pro.p()<4) ? [ele, pro] : null
          }
        }.collectMany{ele,pro->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==321}.findResults{ikp->
            def kp = LorentzVector.withPID(321,*['px','py','pz'].collect{partb.getFloat(it,ikp)})
            return (kp.p()<2.5 && Math.toDegrees(kp.theta())<35) ? [ele,pro,kp] : null
          }
        }.collectMany{ele,pro,kp->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==-321}.findResults{ikm->
            def km = LorentzVector.withPID(-321,*['px','py','pz'].collect{partb.getFloat(it,ikm)})
            return (km.p()<2.5 && Math.toDegrees(km.theta())<40) ? [ele,pro,kp,km] : null
          }
        }.findResults{ele,pro,kp,km->
          def epkpkmX = beam+target-ele-pro-kp-km

          def ekpkmX = beam+target-ele-kp-km
          def epkpX = beam+target-ele-pro-kp
          def epkmX = beam+target-ele-pro-km

          def cplpro = ekpkmX.vect().theta(pro.vect())
          def cplkm = epkpX.vect().theta(km.vect())
          def cplkp = epkmX.vect().theta(kp.vect())
          def pt = Math.sqrt(epkpkmX.px()**2 + epkpkmX.py()**2)
          return (epkpkmX.mass2().abs()<0.06 && cplpro<30 && cplkm<20 && cplkp<20 && pt<0.5) ? [ele,pro,kp,km,epkpkmX,ekpkmX,epkpX,epkmX] : null
        }

      hists.computeIfAbsent('combsize',hcombs).fill(combs.size())
      if(combs.size()==1) {
        def (ele,pro,kp,km,epkpkmX,ekpkmX,epkpX,epkmX) = combs[0]
        def eX = beam+target-ele
        def epX = beam+target-ele-pro
        def vphi = kp+km

        hists.computeIfAbsent('W',hW).fill(eX.mass())
        hists.computeIfAbsent('missemm2',hmissemm2).fill(epkpkmX.e(), epX.mass2())
        hists.computeIfAbsent('phim',hphim).fill(vphi.mass())
        hists.computeIfAbsent('mm2phim',hmissephim).fill(epX.mass2(), vphi.mass())
        hists.computeIfAbsent('missephim',hmissephim).fill(epkpkmX.e(), vphi.mass())
        hists.computeIfAbsent('mismm2phim',hmissmm2phim).fill(epkpkmX.mass2(), vphi.mass())
        hists.computeIfAbsent('missmm2',hmissmm2).fill(epkpkmX.mass2())
      }
    }
  }
}
