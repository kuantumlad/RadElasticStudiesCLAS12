package utils
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.pdg.PDGDatabase


class KinTool{

    //move somplace else

    static def calcQ2(LorentzVector beam, LorentzVector measured_el){
	return -(beam-measured_el).mass2()
    }

    static def calcXb(LorentzVector beam, LorentzVector measured_el){
	return calcQ2(beam,measured_el)/(2.0*PDGDatabase.getParticleMass(2212)*(beam-measured_el).e())
    }

    static def calcT(LorentzVector measured_prot){
	return 2.0*PDGDatabase.getParticleMass(2212)*(measured_prot-PDGDatabase.getParticleMass(2212))
    }

    static def calcNu(LorentzVector beam, LorentzVector measured_el){
	return calcQ2(beam,measured_el)/(2*PDGDatabase.getParticleMass(2212)*calcXb(beam,measured_el))
    }
    
    static def calcY(LorentzVector beam, LorentzVector measured_el){
	return calcNu(beam,measured_el)/beam.e()
    }

    static def calcPhiTrento(LorentzVector beam, LorentzVector measured_el,LorentzVector measured_prot){
	def v3l = beam.vect().cross(ele.vect())
	def v3h = pro.vect().cross((beam-ele).vect())
	def trento = Math.toDegrees(Math.acos(v3l.dot(v3h)/(v3l.mag()*v3h.mag())))

	if(v3l.dot(pro.vect())<0){trento=360-trento}
	return trento
    }

    static def delta_meas_energy( Double beam,  LorentzVector measured_el ){
	def calc_e = beam/(1+ (beam/PDGDatabase.getParticleMass(2212))*(1-Math.cos(measured_el.theta())) );
	return (calc_e  - measured_el.e() )
    }

    static def delta_meas_theta( Double beam,  LorentzVector measured_el ){
	def calc_theta = Math.toDegrees(Math.acos( 1 + (PDGDatabase.getParticleMass(2212)/beam)*( 1 - beam/measured_el.e()) ))
	return Math.toDegrees((calc_theta  - measured_el.theta() ))
    }

    
    

}


