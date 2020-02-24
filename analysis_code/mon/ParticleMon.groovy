package mon
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import groovy.sql.Sql

class ParticleMon {
  def plots = [:].withDefault{[:].withDefault{new H1F(it, "energy", 200,0,10)}}
  def data = [:].withDefault{[:].withDefault{0}}
  def dets = ['FT', 'FD', 'CD']
  def qname = [(-1):'neg', 0:'neut', 1:'pos']
  def pname = [2212:'proton', 11:'electron', 211:'pip', (-211):'pim']

  def processEvent(event) {
    def secs = [:]
    if (event.hasBank("REC::Calorimeter")) {
      def calbank = event.getBank("REC::Calorimeter")
      secs << [calbank.getShort('pindex'),calbank.getByte('detector'),calbank.getByte('sector')].transpose()
        .findAll{it[1].toInteger() == DetectorType.ECAL.getDetectorId()}
        .collectEntries{ [(it[0].toInteger()): it[2]] }
    }
    if (event.hasBank("REC::Scintillator")) {
      def scibank = event.getBank("REC::Scintillator")
      secs << [scibank.getShort('pindex'),scibank.getByte('detector'),scibank.getByte('sector')].transpose()
        .findAll{it[1].toInteger() == DetectorType.FTOF.getDetectorId()}
        .collectEntries{ [(it[0].toInteger()): it[2]] }
    }

    if (event.hasBank("REC::Particle")) {
      def pbank = event.getBank("REC::Particle")

      (0..<pbank.rows()).each{ind ->
        int iq = pbank.getByte('charge', ind)
        int pid = pbank.getInt('pid', ind)
        int idet = pbank.getShort('status',ind).abs()/2000
        /**/

        if(qname[iq]!=null) {
          def name0 = qname[iq] + "@" + dets[idet]
          def name1 = (idet==1) ? ('sec'+secs[ind]) : name0
          if(!name1.contains('null')) {
            data[name0][name1]++
            plots[name0][name1].fill(Math.sqrt(['px','py','pz'].collect{pbank.getFloat(it,ind)}.collect{it*it}.sum()))
          }
        }
        if(pname[pid]!=null) {
          def name0 = pname[pid] + "@" + dets[idet]
          def name1 = (idet==1) ? ('sec'+secs[ind]) : name0
          if(!name1.contains('null')) {
            data[name0][name1]++
            plots[name0][name1].fill(Math.sqrt(['px','py','pz'].collect{pbank.getFloat(it,ind)}.collect{it*it}.sum()))
          }
        }
      }
//    def proton = new Vector3(*['px', 'py', 'pz'].collect{pbank.getFloat(it, proind)})
    }
  }
}

