package mon
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import groovy.sql.Sql
import pid.electron.Electron

class GoodElectron {
  def plots = [:].withDefault{[:].withDefault{new H1F("heth_$it", "electron theta", 90,0,45)}}
  def data = [:].withDefault{[:].withDefault{0}}

  def processEvent(event) {
    if (event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter")) {
      def pbank = event.getBank("REC::Particle")
      def calbank = event.getBank("REC::Calorimeter")

      int eleind = (0..<pbank.rows()).find{pbank.getInt("pid",it)==11 && pbank.getShort("status",it)<0}
      def elesec = (0..<calbank.rows()).collect{
        (calbank.getShort('pindex',it).toInteger() == eleind &&
        calbank.getByte('detector',it).toInteger() == DetectorType.ECAL.getDetectorId()) ? calbank.getByte('sector',it) : null
      }.find()

      if(eleind!=null && elesec!=null) {
        def electron = new Particle(11, *['px', 'py', 'pz'].collect{pbank.getFloat(it, eleind)})
        plots.trigger[elesec].fill(electron.theta())
        data.trigger["sec$elesec"]++
      }
    }
  }
}

