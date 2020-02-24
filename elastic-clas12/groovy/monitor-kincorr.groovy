#!/opt/coatjava/6.3.1/bin/run-groovy
import event.Event
import event.EventConverter
import groovyx.gpars.GParsPool
import java.util.concurrent.ConcurrentHashMap
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.io.hipo.HipoDataSource

import org.jlab.jroot.ROOTFile
import org.jlab.jroot.TNtuple

// Additional class members to make the object
// more useful for CLAS12 analysis.
Particle.metaClass.pindex = null
Particle.metaClass.sphi = null
Particle.metaClass.sector = null

def beam = new Particle(11, 0.0, 0.0, 10.604)
def target = new Particle(2212, 0.0, 0.0, 0.0)

def data_mm2_cut_scaler = 3

cuts = [
    w: [0.8, 1.15],
    high_w:[1.15,999.9],
    w_loose: [0.8, 1.30],
    angle: [178, 185],
    missing_pt: [0.0, 0.2],
    theta_gamma: [0, 3],
    p_ele:[1.5, 10.646],
    missing_mass:[-0.4, 0.4],
    missing_mass_ftof:[-0.1,0.1],
    missing_mass_rad:[-0.01*data_mm2_cut_scaler, 0.01*data_mm2_cut_scaler],
    missing_mass_rad_ftof:[-0.01*data_mm2_cut_scaler ,0.01*data_mm2_cut_scaler]
]

tighter_kin_bounds = [
        theta_ele : [5, 45],
        theta_pro : [5, 90],
        theta_sum : [0, 120],
        p_ele     : [0.1, 10.5],
        p_pro     : [0.1, 5.5],
        w         : [0.6, 4.7],
        x         : [0.0, 1.0],
        phi       : [-30, 330],
        phi_local : [-35,35],
        dp_ele    : [-1, 1],        
        dp_pro    : [-1, 1],
        dp_ele_small    : [-0.2, 0.2],
        dp_pro_small    : [-0.2, 0.2],
        dtheta_ele: [-5, 5],
        dtheta_pro: [-5, 5],
        angle_ep  : [120, 185],
        q2        : [1.2, 4.5],
    missing_pt    : [0, 1],
    e_gamma: [0, 11],
    theta_gamma:[0, 35],
    theta_egamma:[0, 35],
    chi2:[0,10],
    missing_mass_small: [-0.3, 0.3],
    missing_mass: [-1, 1],
    de_beam   : [-2, 2]
   
]

lim = tighter_kin_bounds

theta_bin_ranges = [
    theta_range_ctof_isr : 3,
    tb_ctof_isr : 12,
    theta_range_ftof_isr : 4,
    tb_ftof_isr : 7,
    theta_range_ctof_fsr : 3,
    tb_ctof_fsr : 12,
    theta_range_ftof_fsr : 4,
    tb_ftof_fsr : 7,
    theta_range_ctof_fsr_pro : 5,
    tb_ctof_fsr_pro : 5,
    theta_range_ftof_fsr_pro : 4,
    tb_ftof_fsr_pro : 6
]


def limited_h1 = { title, nbins, lims ->
    new H1F("$title", "$title", nbins, lims[0], lims[1])
}

def limited_h2 = { title, nxbins, nybins, xlims, ylims ->
    new H2F("$title", "$title", nxbins, xlims[0], xlims[1], nybins, ylims[0], ylims[1])
}

def h_theta_bin_elastic = new H1F("elastic_thetabins","elastic_thetabins", 6, 4.0, 16.0)
def h_theta_bin_isr = new H1F("isr_thetabins","isr_thetabins", 5, 7.0, 35.0)
def h_theta_bin_isr_pro = new H1F("isr_thetabins_pro","isr_thetabins_pro", 5, 20.0, 35.0) // same range for fsr
def h_theta_bin_fsr = new H1F("fsr_thetabins","fsr_thetabins", 5, 5.0, 25.0)
def h_phi_bin = new H1F("phibins","phibins", 12, -30.0, 30.0)


histos = new ConcurrentHashMap()
histoBuilders = [
        w        : { title -> limited_h1(title, 200, lim.w) },
        p_ele    : { title -> limited_h1(title, 200, lim.p_ele) },
        p_pro    : { title -> limited_h1(title, 200, lim.p_pro) },
        dp_ele    : { title -> limited_h1(title, 200, lim.dp_ele_small) },
        dp_pro    : { title -> limited_h1(title, 200, lim.dp_pro_small) },
        de_beam  : { title -> limited_h1(title, 200, lim.de_beam) },
        angle_ep : { title -> limited_h1(title, 200, lim.angle_ep) },
        theta_p  : { title -> limited_h1(title, 200, lim.theta_pro) },
        theta_ele: { title -> limited_h1(title, 200, lim.theta_ele) },
        theta_pro: { title -> limited_h1(title, 200, lim.theta_pro) },
        e_gamma: { title -> limited_h1(title, 200, lim.e_gamma) },
        theta_gamma: { title -> limited_h1(title, 200, lim.theta_gamma) },
        missing_mass: { title -> limited_h1(title, 200, lim.missing_mass) },
        missing_mass_small: { title -> limited_h1(title, 200, lim.missing_mass_small) }
]

histoBuilders2 = [
        p_pro_dp         : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.dp_pro) },
        p_ele_dp         : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.dp_ele) },
        p_ele_dp_small   : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.dp_ele_small) },
        phi_ele_dp       : { title -> limited_h2(title, 100, 100, lim.phi_local, lim.dp_ele_small) },
        p_ele_theta      : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.theta_ele) },
        p_pro_theta      : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.theta_pro) },
        theta_ele_dtheta : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.dtheta_ele) },    
        theta_pro_dtheta : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.dtheta_pro) },
        theta_ele_dp     : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.dp_ele) },
        theta_pro_dp     : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.dp_pro) },
        p_ele_dtheta     : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.dtheta_ele) },
        p_pro_dtheta     : { title -> limited_h2(title, 200, 200, lim.p_pro, lim.dtheta_pro) },
        w_q2             : { title -> limited_h2(title, 200, 200, lim.w, lim.q2) },
        x_q2             : { title -> limited_h2(title, 200, 200, lim.x, lim.q2) },
    theta_theta: { title -> limited_h2(title, 200, 200, lim.theta_egamma, lim.theta_gamma) },
    w_theta_sum : { title -> limited_h2(title, 200, 200, lim.w, lim.theta_sum) },
    w_theta_gamma : { title -> limited_h2(title, 200, 200, lim.w, lim.theta_gamma) },
    de_beam_de_beam  : { title -> limited_h2(title, 200, 200, lim.de_beam, lim.de_beam) },
    theta_ele_de_beam: { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.de_beam) },
    theta_pro_de_beam: { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.de_beam) },
    chi2_dp_ele      : { title -> limited_h2(title, 100, 100, lim.chi2, lim.dp_ele) },
    chi2_dp_pro      : { title -> limited_h2(title, 100, 100, lim.chi2, lim.dp_pro) },
]


def loadMomentumCorrectionFiles(path, ftof_or_ctof, sector, n_thetabins){
    println(' Loading Momentum Correction Files ')
    println(' path ' + path + ' detector ' + ftof_or_ctof + ' S ' + sector)
    kin_corr=[:]
    for( nn in 0..n_thetabins){
	kin_corr_constants = []
	def file = new File(path+"p_ele_dp_ele_from_angles_${ftof_or_ctof}_${sector}_thetabin${nn}_isr.txt")
	if( file.exists() ){
	    file.eachLine {  
		line -> println "line : $line"; 
		def arr = line.tokenize(' ')
		println("${arr[0]} and ${arr[1]}")
		kin_corr_constants.add(Double.valueOf(arr[1]))
	    }
	}
	kin_corr.put(nn,kin_corr_constants)	
    }
    return kin_corr
}

def momentumCorrection(delta_p, p_corr_parameters){
    //println(p_corr_parameters)
    def a = p_corr_parameters[0]
    def b = p_corr_parameters[1]
    def c = p_corr_parameters[2]    
    return c*delta_p*delta_p + b*delta_p + a
}

def getGeneratedSector(phi){
    return Math.ceil(phi / 60) + 3
}

def shiftPhi(phi) {
    return (phi > 150) ? phi - 180 : phi + 180
}

def relativePhi(phi, sector) {
    if (sector == 4) {
        return phi
    } else if (sector == 5) {
        return phi - 60
    } else if (sector == 6) {
        return phi - 120
    } else if (sector == 1) {
        return phi - 180
    } else if (sector == 2) {
        return phi - 240
    } else if (sector == 3) {
        return phi - 300
    }
}

def getKin(beam, target, electron) {

    def missing = new Particle(beam)
    missing.combine(target, 1)
    missing.combine(electron, -1)
    def w = missing.mass()

    def q = new Particle(beam)
    q.combine(electron, -1)
    def q2 = -1 * q.mass2()

    def nu = beam.e() - electron.e()
    def y = nu / beam.e()
    def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))

    return [x: x, y: y, w: w, nu: nu, q2: q2]
}

def getPKin(beam, target, electron, proton) {

    def missing = new Particle(beam)
    missing.combine(target, 1)
    missing.combine(electron, -1)
    def w = missing.mass()
    missing.combine(proton, -1)
    def missing_mass = missing.mass2()

    def q = new Particle(beam)
    q.combine(electron, -1)
    def q2 = -1 * q.mass2()

    def nu = beam.e() - electron.e()
    def y = nu / beam.e()
    def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))

    def zaxis = new Vector3(0, 0, 1)
    def enorm = electron.vector().vect().cross(zaxis)
    def pnorm = proton.vector().vect().cross(zaxis)
    def phi = enorm.theta(pnorm)
    //def phi = Math.toDegrees(electron.phi() - proton.phi())
    def missing_pt = Math.sqrt(missing.px()**2 + missing.py()**2)
    def theta_gamma = Math.toDegrees(missing.theta())
    def norm = missing.vector().vect().mag() * electron.vector().vect().mag()
    def theta_egamma = Math.toDegrees(Math.acos(missing.vector().vect().dot(electron.vector().vect()) / norm))

    return [x: x, y: y, w: w, nu: nu, q2: q2, angle: phi,
            missing_mass: missing_mass, missing_energy: missing.e(),
	    missing_pt: missing_pt, theta_gamma:theta_gamma, 
	    theta_egamma:theta_egamma]
}


//use this function for predicting the FSR proton P using Beam and Theta_pr
def predictProtonMathematica(beam,theta){
    //doesn't work
    //println("Calculating predicted proton momentum using proton angle")
    def M = PDGDatabase.getParticleMass(2212)
    // get the momentum of the proton
    def num = 4*beam.e()*M*(beam.e() + M)*Math.cos(theta)
    def den = beam.e()**2 + 4*beam.e()*M + 2*(M**2) - 2*(beam.e()**2)*Math.cos(2*theta)
    def pred_pr_mntm = num/den
    //println(' pred mntm ' + pred_pr_mntm ) 
   return pred_pr_mntm // predicted proton momentum using beam energy and proton theta 
}

//predict proton angle using electron angle
def predictProtonAngleMathematics(beam, theta){
    println("predict proton theta using beam and electron angles")
    //get the angle of the proton
    def c0 = (beam.e() + M)/M
    def c1 = Math.tan(theta/2)
    def a0 = (Math.PI/2.0) - Math.atan(co*c1)
    def pred_pr_theta = Math.toDegrees(a0)
 
    println(' proton pred angle ' + pred_pr_theta)
    return pred_pr_theta
}

//this works
def predictElectronMathematica(pro){
    // This one actually works! 
    def M = PDGDatabase.getParticleMass(2212)
    def Pp = pro.p()
    def Ep = Math.sqrt(M**2 + Pp**2)
    def beta = pro.theta()
    def cbeta = Math.cos(beta)
    def c2beta = Math.cos(2 * beta)
    def c3beta = Math.cos(3 * beta)
    def sbeta = Math.sin(-beta)
    
    // Angle 
    def expr_num = 2 * M**2 - 2 * M * Ep - Pp * (M + Ep) * cbeta
    expr_num += 2 * M * (M + Ep) * c2beta + M * Pp * c3beta + Pp * Ep * c3beta
    def expr_den = 2 * (2 * M**2 + Pp**2 - Pp**2 * c2beta)
    def alpha = Math.acos(-1 * expr_num/expr_den)
    
    // Mom
    expr_num = -2 * M * (-M + Ep)  * cbeta + Pp * (M + Ep) + (M - Ep) * c2beta
    expr_den = 4 * M * cbeta - 2 * Pp * sbeta**2
    def kprime = expr_num / expr_den

    return [momentum:kprime, theta:alpha]
}

//use for predicing Ps in the initial state radiation case (ISR)
//requires only angles
//this works
def predictElectronMomentum(ele, pro){
    def M = PDGDatabase.getParticleMass(2212)
    return M * Math.cos(ele.theta()/2 + pro.theta()) / Math.sin(ele.theta()/2) / Math.sin(ele.theta() + pro.theta())
}

//this works
def predictProtonMomentum(ele, pro){
    def M = PDGDatabase.getParticleMass(2212)
    return 2 * M * Math.cos(ele.theta()/2) * Math.cos(ele.theta()/2 + pro.theta()) / Math.sin(pro.theta()) / Math.sin(ele.theta() + pro.theta())    
}

//use for FSR case when finding  the deltas in P and Theta
// this works!
def predictElasticBasedOnElectronAngle(beam, alpha) {
    // alpha - electron angle in radians (known)
    // beta  - proton angle in radians (inferred)
    //
    // Don't confuse the beta in this function with the
    // physics beta, v/c.

    // This is stupid, but it makes the next lines shorter
    def proton_mass = PDGDatabase.getParticleMass(2212)
    def cos_alpha = Math.cos(alpha)
    def sin_alpha = Math.sin(alpha)

    // Kinematics prediction
    def pred_ele_p = beam.e() / (1 + (beam.e() / proton_mass) * (1 - cos_alpha))
    def pred_pro_p = Math.sqrt(beam.e()**2 - 2 * beam.e() * pred_ele_p * cos_alpha + pred_ele_p**2)
    def pred_beta = Math.asin((pred_ele_p * sin_alpha) / pred_pro_p)

    return [pred_ele_p, pred_beta, pred_pro_p]
}

//this function checks out - it is okay
def predictBeamEnergy(ele, pro){
    def pred_e_beam = ele.p() / (1 + (ele.p() / PDGDatabase.getParticleMass(2212)) * (Math.cos(ele.theta()) - 1))
    def a0 = 1 - 1 / (Math.cos(ele.theta()) - Math.sin(ele.theta()) / Math.tan(-1 * pro.theta()))
    def a1 = Math.sin(ele.theta()) / Math.sin(ele.theta() + pro.theta())
    def pred_e_beam_from_angles = 2 * PDGDatabase.getParticleMass(2212) * a0 / (a1**2 - a0**2)
    return [pred_e_beam, pred_e_beam_from_angles]
}

//////////////////////////////////////////////////////
// load momentum correction factors
def path_to_corr_files = '/w/hallb-scifs17exp/clas12/bclary/CLAS12/david_elastic/elastic-clas12/python/fit_parameters/'
def el_s1_prftof = loadMomentumCorrectionFiles(path_to_corr_files, 'FTOF', '1', 7)
def el_s2_prftof = loadMomentumCorrectionFiles(path_to_corr_files, 'FTOF', '2', 7)
def el_s3_prftof = loadMomentumCorrectionFiles(path_to_corr_files, 'FTOF', '3', 7)
def el_s4_prftof = loadMomentumCorrectionFiles(path_to_corr_files, 'FTOF', '4', 7)
def el_s5_prftof = loadMomentumCorrectionFiles(path_to_corr_files, 'FTOF', '5', 7)
def el_s6_prftof = loadMomentumCorrectionFiles(path_to_corr_files, 'FTOF', '6', 7)

def el_s1_pr_ctof = loadMomentumCorrectionFiles(path_to_corr_files, 'CTOF', '1', 7)
def el_s2_pr_ctof = loadMomentumCorrectionFiles(path_to_corr_files, 'CTOF', '2', 7)
def el_s3_pr_ctof = loadMomentumCorrectionFiles(path_to_corr_files, 'CTOF', '3', 7)
def el_s4_pr_ctof = loadMomentumCorrectionFiles(path_to_corr_files, 'CTOF', '4', 7)
def el_s5_pr_ctof = loadMomentumCorrectionFiles(path_to_corr_files, 'CTOF', '5', 7)
def el_s6_pr_ctof = loadMomentumCorrectionFiles(path_to_corr_files, 'CTOF', '6', 7)

def el_mntm_corr_par = [el_s1_prftof,
			el_s2_prftof,
			el_s3_prftof,
			el_s4_prftof,
			el_s5_prftof,
			el_s6_prftof
                       ] 

def el_mntm_corr_par_pr_ctof = [el_s1_pr_ctof,
				el_s2_pr_ctof,
				el_s3_pr_ctof,
				el_s4_pr_ctof,
				el_s5_pr_ctof,
				el_s6_pr_ctof
                               ] 


GParsPool.withPool 16, {
    args.eachParallel { filename ->


	def reader = new HipoDataSource()
        reader.open(filename)

        def eventIndex = 0
        while (reader.hasEvent()) {
            if (eventIndex % 5000 == 0) {
                println("Processing " + eventIndex)
            }
	    //if( eventIndex == 2000000 ) break

            def dataEvent = reader.getNextEvent()
            def event = EventConverter.convert(dataEvent)

	    // Reconstructed (for data and simulation)
            (0..<event.npart).find {
                event.pid[it] == 11 && event.status[it] < 0 && event.p[it] > cuts.p_ele[0]
            }?.each { idx ->
                def sector = event.dc_sector[idx]
                def ele = new Particle(11, event.px[idx], event.py[idx], event.pz[idx])
                ele.sector = sector
                ele.pindex = idx
                ele.sphi = shiftPhi(Math.toDegrees(ele.phi()))

                def kin = getKin(beam, target, ele)

		// Believe the electron angular measurement and infer the rest of the
                // event kinematics based on the assumption that it's an elastic scattering.
                def (pred_ele_p_ang, pred_pro_theta_ang, pred_pro_p_ang) = predictElasticBasedOnElectronAngle(beam, ele.theta())
		             
                (0..<event.npart).findAll { event.pid[it] == 2212 }.each {
		    def sector_pro = event.dc_sector[it]
                    def pro = new Particle(2212, event.px[it], event.py[it], event.pz[it])
                    pro.pindex = it
                    pro.sphi = shiftPhi(Math.toDegrees(pro.phi()))
                    def pkin = getPKin(beam, target, ele, pro)

		    
		    def ctof = event.ctof_status.contains(it).findResult{ stat -> stat ? "CTOF" : "FTOF"}


		    //make predictions here for Resolutions 

		    def (pred_e_beam, pred_e_beam_from_angles) = predictBeamEnergy(ele, pro)		  
		    def effective_beam = new Particle(11, 0.0, 0.0, pred_e_beam_from_angles)
                    def effective_pkin = getPKin(effective_beam, target, ele, pro)

		    //trust pro P and Theta and Beam - use for FSR
		    def pred_ele = predictElectronMathematica(pro)
		    //predict proton momentum using beam and theta of proton
		    //probably not okay
		    //def pred_pro = predictProton(ele)
		    def pred_pro_p_beam_prang = predictProtonMathematica(beam, Math.toDegrees(pro.theta()))
		    
		    //momentum calculated based on scattering angles - use for ISR 
		    def pred_ele_p = predictElectronMomentum(ele, pro)
		    def pred_pro_p = predictProtonMomentum(ele, pro)


		    def pass_high_w = pkin.w > cuts.high_w[0]
		    def pass_angle_ep = pkin.angle > cuts.angle[0] && pkin.angle < cuts.angle[1]
		    def pass_w_elastic = pkin.w < cuts.high_w[0]
		    def pass_theta_gamma = pkin.theta_gamma < cuts.theta_gamma[1] // IRS CUT
		    def pass_theta_egamma = pkin.theta_egamma < cuts.theta_gamma[1] // FRS CUT
		    
		    
		    def pass_missing_mass_ftof = pkin.missing_mass > cuts.missing_mass_ftof[0] && pkin.missing_mass < cuts.missing_mass_ftof[1]
		    def pass_missing_mass = pkin.missing_mass > cuts.missing_mass[0] && pkin.missing_mass < cuts.missing_mass[1] 
		    pass_missing_mass = (pass_missing_mass && ctof == "CTOF") || (pass_missing_mass_ftof && ctof == "FTOF")

		    def isr_fsr_flag = null		    
  		    //decide which set event belongs to - initial state radiation or final state radiation
		    if( pass_theta_gamma ){
			isr_fsr_flag = 'isr'			
		    }
		    else if( pass_theta_egamma ){
			isr_fsr_flag = 'fsr'
		    }
		    
		    // define new missing mass based on isr or fsr
		    def cor_mm2 = 0
		    def pass_missing_mass_rad = false
		    if( isr_fsr_flag == 'isr' ){
			cor_mm2 = effective_pkin.missing_mass
 			def pass_missing_mass_rad_ftof = cor_mm2 > cuts.missing_mass_rad_ftof[0] && cor_mm2 < cuts.missing_mass_rad_ftof[1]
			pass_missing_mass_rad = cor_mm2 > cuts.missing_mass_rad[0] && cor_mm2 < cuts.missing_mass_rad[1]
			pass_missing_mass_rad = (pass_missing_mass_rad && ctof == "CTOF") || (pass_missing_mass_rad_ftof && ctof == "FTOF" )
		    }
		    else{
			cor_mm2 = pkin.missing_mass
			pass_missing_mass_rad = true
		    }




		    // raw distributions
                    histos.computeIfAbsent('w_' + ctof + '_' + sector, histoBuilders.w).fill(pkin.w)
                    histos.computeIfAbsent('angle_ep_' + ctof, histoBuilders.angle_ep).fill(pkin.angle)
                    histos.computeIfAbsent('theta_gamma_' + ctof, histoBuilders.theta_gamma).fill(pkin.theta_gamma)		   
                    histos.computeIfAbsent('theta_egamma_' + ctof, histoBuilders.theta_gamma).fill(pkin.theta_egamma)
		    histos.computeIfAbsent('theta_egamma_theta_gamma_' + ctof, histoBuilders2.theta_theta).fill(pkin.theta_gamma, pkin.theta_egamma)		    
		    histos.computeIfAbsent('missing_mass_' + ctof, histoBuilders.missing_mass).fill(pkin.missing_mass)
		    histos.computeIfAbsent('w_theta_sum_' + ctof , histoBuilders2.w_theta_sum).fill(
			pkin.w, Math.toDegrees(ele.theta() + pro.theta())
		    )
		    histos.computeIfAbsent('p_ele_theta_ele_' + ctof, histoBuilders2.p_ele_theta).fill(
			ele.p(), Math.toDegrees(ele.theta()))
		    histos.computeIfAbsent('p_pro_theta_pro_' + ctof, histoBuilders2.p_pro_theta).fill(
			pro.p(), Math.toDegrees(pro.theta()))
		    histos.computeIfAbsent('w_theta_sum_' + ctof + '_' + isr_fsr_flag , histoBuilders2.w_theta_sum).fill(
			pkin.w, Math.toDegrees(ele.theta() + pro.theta())
		    )
		    		    		    

		    //println(' >>  ' + isr_fsr_flag + ' missing mass ' + pkin.missing_mass + ' cor_mm2 ' + cor_mm2 )
		    if( isr_fsr_flag != null ) histos.computeIfAbsent('effective_missing_mass_' + ctof, histoBuilders.missing_mass_small).fill(cor_mm2)
		    
		    if( pass_w_elastic ){
			histos.computeIfAbsent('angle_ep_pass_w_elastic_' + ctof, histoBuilders.angle_ep).fill(pkin.angle)
			histos.computeIfAbsent('w_theta_sum_pass_w_elastic_' + ctof , histoBuilders2.w_theta_sum).fill(
			    pkin.w, Math.toDegrees(ele.theta() + pro.theta())
			)
			histos.computeIfAbsent('theta_gamma_pass_w_elastic_' + ctof + '_' + sector, histoBuilders.theta_gamma).fill(pkin.theta_gamma)
  			histos.computeIfAbsent('theta_egamma_pass_w_elastic_' + ctof + '_' + sector, histoBuilders.theta_gamma).fill(pkin.theta_egamma)
		    }
 		    if( pass_w_elastic && pass_missing_mass && pass_angle_ep){
			histos.computeIfAbsent('w_pass_all_' + ctof + '_' + sector + '_pass_all_elastic', histoBuilders.w).fill(pkin.w)
			histos.computeIfAbsent('p_ele_theta_ele_' + ctof + '_' + sector + '_pass_all_elastic', histoBuilders2.p_ele_theta).fill(
			    ele.p(), Math.toDegrees(ele.theta()))
			histos.computeIfAbsent('p_pro_theta_pro_' + ctof + '_' + sector_pro + '_pass_all_elastic', histoBuilders2.p_pro_theta).fill(
			    pro.p(), Math.toDegrees(pro.theta()))
			histos.computeIfAbsent('w_theta_sum_pass_all_' + ctof + '_pass_all_elastic', histoBuilders2.w_theta_sum).fill(
			    pkin.w, Math.toDegrees(ele.theta() + pro.theta()))
			    			    									
 			histos.computeIfAbsent('p_ele_' + ctof + '_' + sector + '_pass_all_elastic', histoBuilders.p_ele).fill(ele.p())
			histos.computeIfAbsent('theta_ele_' + ctof + '_' + sector + '_pass_all_elastic', histoBuilders.theta_ele).fill(Math.toDegrees(ele.theta()))
			histos.computeIfAbsent('p_pro_' + ctof + '_' + sector_pro + '_pass_all_elastic', histoBuilders.p_pro).fill(pro.p())
			histos.computeIfAbsent('theta_pro_' + ctof + '_' + sector_pro + '_pass_all_elastic', histoBuilders.theta_pro).fill(Math.toDegrees(pro.theta()))

			//okay
			histos.computeIfAbsent('theta_ele_de_beam_from_angles_' + ctof + '_' + ele.sector + '_elastic', histoBuilders2.theta_ele_de_beam).fill(
			    Math.toDegrees(ele.theta()), beam.e() - pred_e_beam_from_angles)
			//okay
			histos.computeIfAbsent('theta_pro_de_beam_from_angles_' + ctof + '_' + sector_pro + '_elastic', histoBuilders2.theta_pro_de_beam).fill(
			    Math.toDegrees(pro.theta()), beam.e() - pred_e_beam_from_angles)			

			//okay
			histos.computeIfAbsent('p_ele_dp_ele_from_angles_' + ctof + '_' + sector + '_pass_all_elastic', histoBuilders2.p_ele_dp).fill(
 			    ele.p(), pred_ele_p - ele.p())
			//okay
			histos.computeIfAbsent('theta_ele_dtheta_ele_from_beam_pangle_' + ctof + '_' + sector + '_pass_all_elastic', histoBuilders2.theta_ele_dtheta).fill(
			    Math.toDegrees(ele.theta()), Math.toDegrees(pred_ele.theta - ele.theta()))

			//okay
 			histos.computeIfAbsent('p_pro_dp_pro_from_beam_eangle_' + ctof + '_' + sector_pro + '_pass_all_elastic', histoBuilders2.p_pro_dp).fill(
			    pro.p(), (pred_pro_p_ang - pro.p()))
			
			//okay
			histos.computeIfAbsent('theta_pro_dtheta_pro_from_beam_eangle_' + ctof + '_' + sector_pro + '_pass_all_elastic', histoBuilders2.theta_pro_dtheta).fill(
			    Math.toDegrees(pro.theta()), Math.toDegrees(pred_pro_theta_ang - pro.theta()))

			//okay
			histos.computeIfAbsent('p_ele_dtheta_ele_from_angles_' + ctof + '_' + sector + '_elastic', histoBuilders2.p_ele_dtheta).fill(
 			    ele.p(), Math.toDegrees(pred_ele.theta - ele.theta()))
			//okay
			histos.computeIfAbsent('theta_ele_dp_ele_from_angles_' + ctof + '_' + sector + '_elastic', histoBuilders2.theta_ele_dp).fill(
 			    Math.toDegrees(ele.theta()), pred_ele_p - ele.p())
			// look at delta p per theta bin

 			def theta_bin_elastic = h_theta_bin_elastic.getXaxis().getBin(Math.toDegrees(ele.theta()))
 			def rel_el_phi = relativePhi(ele.sphi, sector)
			//okay
			histos.computeIfAbsent('p_ele_dp_ele_from_angles_' + ctof + '_' + sector + '_thetabin' + theta_bin_elastic + '_pass_all_elastic', histoBuilders2.p_ele_dp).fill(
 			    ele.p(), pred_ele_p - ele.p())
			histos.computeIfAbsent('phi_local_ele_dpp_ele_from_angles_' + ctof + '_' + sector + '_thetabin' + theta_bin_elastic + '_pass_all_elastic', histoBuilders2.phi_ele_dp).fill(
			    rel_el_phi, (pred_ele_p - ele.p())/ele.p())
			histos.computeIfAbsent('phi_local_ele_dp_ele_from_angles_' + ctof + '_' + sector + '_thetabin' + theta_bin_elastic + '_pass_all_elastic', histoBuilders2.phi_ele_dp).fill(
			    rel_el_phi, (pred_ele_p - ele.p()))
		    }
		

    		    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		    /////////////////////////////////     Second part of analysis is to look at radiated elastic    ///////////////////////////////////////////////////////////////////
		    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		    
		    if( pass_angle_ep ){ // want to see it with just angle ep cut, no isr or fsr cuts
			histos.computeIfAbsent('w_theta_sum_' + ctof + '_pass_angle_ep', histoBuilders2.w_theta_sum).fill(
			    pkin.w, Math.toDegrees(ele.theta() + pro.theta()))
		    }
		    if( pass_angle_ep && pass_high_w ){ // want to see it with just angle ep cut ang high W cut
			histos.computeIfAbsent('w_theta_sum_' + ctof + '_pass_angle_ep_highw', histoBuilders2.w_theta_sum).fill(
			    pkin.w, Math.toDegrees(ele.theta() + pro.theta()))
		    }
		    if( pass_angle_ep && (pass_missing_mass_rad || pass_missing_mass )){ // show before ISR or FSR cuts
			histos.computeIfAbsent('w_theta_sum_' + ctof + '_pass_angle_ep_mm2_rad_elastic', histoBuilders2.w_theta_sum).fill(
			    pkin.w, Math.toDegrees(ele.theta() + pro.theta()))
		    }

		    
		    if (pass_theta_gamma ){
			histos.computeIfAbsent('w_pass_theta_gamma_' + ctof + '_' + sector, histoBuilders.w).fill(pkin.w)
			histos.computeIfAbsent('angle_ep_pass_theta_gamma_' + ctof + '_' + sector, histoBuilders.angle_ep).fill(pkin.angle)
		    }

		    if (pass_theta_egamma ){
			histos.computeIfAbsent('w_pass_theta_egamma_' + ctof + '_' + sector, histoBuilders.w).fill(pkin.w)
			histos.computeIfAbsent('angle_ep_pass_theta_egamma_' + ctof + '_' + sector, histoBuilders.angle_ep).fill(pkin.angle)
		    }

		    //plane angle cut
		    if (pkin.angle > cuts.angle[0] && pkin.angle < cuts.angle[1]){
			histos.computeIfAbsent('w_pass_angle_' + ctof + '_' + sector, histoBuilders.w).fill(pkin.w)
 			histos.computeIfAbsent('w_pass_angle_'  + ctof, histoBuilders.w).fill(pkin.w)
			histos.computeIfAbsent('theta_gamma_pass_angle_' + ctof + '_' + sector, histoBuilders.theta_gamma).fill(pkin.theta_gamma)
  			histos.computeIfAbsent('theta_egamma_pass_angle_' + ctof + '_' + sector, histoBuilders.theta_gamma).fill(pkin.theta_egamma)
			histos.computeIfAbsent('theta_e_theta_gamma_pass_angle_' + ctof + '_' + sector, histoBuilders2.theta_theta).fill(
			    Math.toDegrees(ele.theta()), pkin.theta_gamma
			)
			histos.computeIfAbsent('theta_e_theta_egamma_pass_angle_' + ctof + '_' + sector, histoBuilders2.theta_theta).fill(
			    Math.toDegrees(ele.theta()), pkin.theta_egamma
			)
			
			//def pass_missing_mass = pkin.missing_mass > cuts.missing_mass[0] && pkin.missing_mass < cuts.missing_mass[1]
			if( isr_fsr_flag != null && pass_high_w && pass_angle_ep ){
			    histos.computeIfAbsent('w_theta_sum_' + ctof + '_rad_elastic1', histoBuilders2.w_theta_sum).fill(
				pkin.w, Math.toDegrees(ele.theta() + pro.theta()))
			    if( pass_missing_mass_rad ){
				histos.computeIfAbsent('w_theta_sum_' + ctof + '_rad_elastic2', histoBuilders2.w_theta_sum).fill(
				    pkin.w, Math.toDegrees(ele.theta() + pro.theta()))
			    }
			}
			

			if( pass_high_w && pass_angle_ep){ //plot gammas with no MM2 otherwise it is indirectly cutting on theta gamma
			    histos.computeIfAbsent('theta_gamma_' + ctof + '_rad_elastic', histoBuilders.theta_gamma).fill(pkin.theta_gamma)
  			    histos.computeIfAbsent('theta_egamma_' + ctof + '_rad_elastic', histoBuilders.theta_gamma).fill(pkin.theta_egamma)			    		
			    histos.computeIfAbsent('theta_gamma_w_' + ctof + '_rad_elastic', histoBuilders2.w_theta_gamma).fill(pkin.w, pkin.theta_gamma)
			    histos.computeIfAbsent('theta_egamma_w_' + ctof + '_rad_elastic', histoBuilders2.w_theta_gamma).fill(pkin.w, pkin.theta_egamma)
			    histos.computeIfAbsent('theta_egamma_theta_gamma_' + ctof + '_rad_elastic', histoBuilders2.theta_theta).fill(pkin.theta_gamma, pkin.theta_egamma)
			}
								
 			// These are ISR events and FSR events			

			if (isr_fsr_flag != null && pass_missing_mass_rad && pass_high_w && pass_angle_ep) {
			    histos.computeIfAbsent('p_ele_theta_ele_' + ctof + '_' + sector + '_rad_elastic', histoBuilders2.p_ele_theta).fill(
				ele.p(), Math.toDegrees(ele.theta()))
			    histos.computeIfAbsent('p_pro_theta_pro_' + ctof + '_' + sector_pro + '_rad_elastic', histoBuilders2.p_pro_theta).fill(
				pro.p(), Math.toDegrees(pro.theta()))

 			    histos.computeIfAbsent('w_pass_all_' + ctof + '_' + sector + '_' + isr_fsr_flag, histoBuilders.w).fill(pkin.w)
			    histos.computeIfAbsent('p_ele_theta_ele_' + ctof + '_' + sector + '_' + isr_fsr_flag, histoBuilders2.p_ele_theta).fill(
			    ele.p(), Math.toDegrees(ele.theta()))
			    histos.computeIfAbsent('p_pro_theta_pro_' + ctof + '_' + sector_pro + '_' + isr_fsr_flag, histoBuilders2.p_pro_theta).fill(
			    pro.p(), Math.toDegrees(pro.theta()))
			    histos.computeIfAbsent('w_theta_sum_pass_all_' + ctof + '_' + isr_fsr_flag, histoBuilders2.w_theta_sum).fill(
				pkin.w, Math.toDegrees(ele.theta() + pro.theta())
			    )


			    histos.computeIfAbsent('p_ele_' + ctof + '_' + sector + '_' + isr_fsr_flag, histoBuilders.p_ele).fill(ele.p())
			    histos.computeIfAbsent('theta_ele_' + ctof + '_' + sector + '_' + isr_fsr_flag, histoBuilders.theta_ele).fill(Math.toDegrees(ele.theta()))
			    histos.computeIfAbsent('p_pro_' + ctof + '_' + sector_pro + '_' + isr_fsr_flag, histoBuilders.p_pro).fill(pro.p())
			    histos.computeIfAbsent('theta_pro_' + ctof + '_' + sector_pro + '_' + isr_fsr_flag, histoBuilders.theta_pro).fill(Math.toDegrees(pro.theta()))



			    //define the theta bin size(theta_range) and number of bins based on El./Proton and ISR and FSR in FTOF/CTOF
			    // determine the theta and phi bin the electron and proton are in

 			    def rel_el_phi = relativePhi(ele.sphi, sector)
			    def rel_pro_phi = relativePhi(pro.sphi, sector_pro)

			    def theta_bin=0
			    def theta_bin_isr = h_theta_bin_isr.getXaxis().getBin(Math.toDegrees(ele.theta()))
			    def theta_bin_fsr = h_theta_bin_fsr.getXaxis().getBin(Math.toDegrees(ele.theta()))
			    
			    if( isr_fsr_flag == 'isr' ) theta_bin = theta_bin_isr
			    else if(isr_fsr_flag == 'fsr') theta_bin = theta_bin_fsr
				
			    //calculate new momentum values
			    def new_p = ele.p()
			    def delta_p = pred_ele.momentum - ele.p()
			    if( isr_fsr_flag == 'isr' && ctof == 'FTOF' && theta_bin_isr >= 0 && theta_bin_isr < 5 ){ 				
				def corr_par = el_mntm_corr_par[sector-1]
				def abc_par = corr_par[theta_bin_isr]
				def delta_p_cor = momentumCorrection(ele.p(), el_mntm_corr_par[sector-1][theta_bin_isr])
				new_p = ele.p()*(1 - delta_p_cor) 
 				//println(' old p ' + ele.p() + ' new p ' + new_p)
				histos.computeIfAbsent('p_ele_dp_ele_from_angles_' + ctof + '_' + sector + '_thetabin' + theta_bin + '_' + isr_fsr_flag + '_corr', histoBuilders2.p_ele_dp_small).fill(
				ele.p(), new_p - ele.p())

			    }
			    if( isr_fsr_flag == 'isr' && ctof == 'CTOF' && theta_bin_fsr >= 0 && theta_bin_fsr < 5 ){
				def corr_par_ctof = el_mntm_corr_par_pr_ctof[sector-1]
				def abc_par = corr_par_ctof[theta_bin_fsr]
				def delta_p_corr = momentumCorrection(ele.p(), abc_par)
				new_p = ele.p() * (1 - delta_p_corr)
				histos.computeIfAbsent('p_ele_dp_ele_from_angles_' + ctof + '_' + sector + '_thetabin' + theta_bin + '_' + isr_fsr_flag + '_corr', histoBuilders2.p_ele_dp_small).fill(
				ele.p(), new_p - ele.p())
			    }
			    
			    
			    def theta_bin_pro = h_theta_bin_isr_pro.getXaxis().getBin(Math.toDegrees(pro.theta()))
			    
			    def phi_bin = h_phi_bin.getXaxis().getBin(Math.toDegrees(ele.phi()))
			    def phi_bin_pro = h_phi_bin.getXaxis().getBin(Math.toDegrees(pro.phi()))
			    
			    

			    // David - only use the _from_angles_ histograms
			    //okay
  			    histos.computeIfAbsent('p_ele_dp_ele_' + ctof + '_' + sector + '_' + isr_fsr_flag, histoBuilders2.p_ele_dp).fill(ele.p(), pred_ele.momentum - ele.p())
			    //okay
			    histos.computeIfAbsent('p_ele_dp_ele_from_angles_' + ctof + '_' + sector + '_' + isr_fsr_flag, histoBuilders2.p_ele_dp).fill(
 				ele.p(), pred_ele_p - ele.p())
			    histos.computeIfAbsent('p_ele_dtheta_ele_from_angles_' + ctof + '_' + sector + '_' + isr_fsr_flag, histoBuilders2.p_ele_dtheta).fill(
 				ele.p(), Math.toDegrees(pred_ele.theta - ele.theta()))
			    histos.computeIfAbsent('theta_ele_dp_ele_from_angles_' + ctof + '_' + sector + '_' + isr_fsr_flag, histoBuilders2.theta_ele_dp).fill(
 				Math.toDegrees(ele.theta()), pred_ele_p - ele.p())
 			    histos.computeIfAbsent('p_ele_dp_ele_from_angles_' + ctof + '_' + sector + '_thetabin' + theta_bin + '_' + isr_fsr_flag, histoBuilders2.p_ele_dp).fill(
				ele.p(), (pred_ele_p - ele.p()))

			    //okay
			    histos.computeIfAbsent('p_ele_dpp_ele_from_angles_' + ctof + '_' + sector + '_' + isr_fsr_flag, histoBuilders2.p_ele_dp).fill(
				ele.p(), (pred_ele_p - ele.p())/ele.p())
 			    // Check delta p/p vs p in different theta bins
			    //okay
 			    histos.computeIfAbsent('p_ele_dpp_ele_from_angles_' + ctof + '_' + sector + '_thetabin' + theta_bin + '_' + isr_fsr_flag, histoBuilders2.p_ele_dp).fill(
				ele.p(), (pred_ele_p - ele.p())/ele.p())
 			    histos.computeIfAbsent('phi_local_ele_dpp_ele_from_angles_' + ctof + '_' + sector + '_thetabin' + theta_bin + '_' + isr_fsr_flag, histoBuilders2.phi_ele_dp).fill(
				rel_el_phi, (pred_ele_p - ele.p())/ele.p())
			    histos.computeIfAbsent('phi_local_ele_dp_ele_from_angles_' + ctof + '_' + sector + '_thetabin' + theta_bin + '_' + isr_fsr_flag, histoBuilders2.phi_ele_dp).fill(
				rel_el_phi, (pred_ele_p - ele.p()))


			    //okay
 			    //histos.computeIfAbsent('p_ele_dpp_ele_from_angles_' + ctof + '_' + sector + '_thetabin' + theta_bin + '_phibin' + phi_bin + '_' + isr_fsr_flag, histoBuilders2.p_ele_dp).fill(ele.p(), (pred_ele_p - ele.p())/ele.p())
			    //okay
			    histos.computeIfAbsent('dp_ele_from_angles_' + ctof + '_' + sector + '_thetabin' + theta_bin + '_' + isr_fsr_flag, histoBuilders.dp_ele).fill((pred_ele_p - ele.p())/ele.p())
			    //okay
			    //histos.computeIfAbsent('dp_ele_from_angles_' + ctof + '_' + sector + '_thetabin' + theta_bin + '_phibin' + phi_bin + '_' + isr_fsr_flag, histoBuilders.dp_ele).fill((pred_ele_p - ele.p())/ele.p())

 			    //okay
			    histos.computeIfAbsent('theta_ele_dtheta_ele_from_beam_pangle_' + ctof + '_' + sector + '_' + isr_fsr_flag, histoBuilders2.theta_ele_dtheta).fill(
				Math.toDegrees(ele.theta()), Math.toDegrees(pred_ele.theta - ele.theta()))
			    //okay
			    histos.computeIfAbsent('theta_ele_dthetath_ele_from_beam_pangle_' + ctof + '_' + sector + '_' + isr_fsr_flag, histoBuilders2.theta_ele_dtheta).fill(
				Math.toDegrees(ele.theta()), Math.toDegrees(pred_ele.theta - ele.theta())/Math.toDegrees(ele.theta()))
			    			   
			    //okay
			    histos.computeIfAbsent('p_pro_dp_pro_' + ctof + '_' + sector_pro + '_' + isr_fsr_flag, histoBuilders2.p_pro_dp).fill(pro.p(), pred_pro_p - pro.p())
			    histos.computeIfAbsent('p_pro_dp_pro_' + ctof + '_' + sector_pro + '_thetabin' + theta_bin_pro + '_' + isr_fsr_flag, histoBuilders2.p_pro_dp).fill(pro.p(), pred_pro_p - pro.p())
			    //okay
			    histos.computeIfAbsent('p_pro_dpp_pro_' + ctof + '_' + sector_pro + '_' + isr_fsr_flag, histoBuilders2.p_pro_dp).fill(pro.p(), (pred_pro_p - pro.p())/pro.p())
			    //okay
			    histos.computeIfAbsent('p_pro_dpp_pro_' + ctof + '_' + sector_pro + '_thetabin' + theta_bin_pro + '_' + isr_fsr_flag, histoBuilders2.p_pro_dp).fill(pro.p(), (pred_pro_p - pro.p())/pro.p())
			    //okay
 			    histos.computeIfAbsent('p_pro_dp_pro_from_beam_eangle_' + ctof + '_' + sector_pro + '_' + isr_fsr_flag, histoBuilders2.p_pro_dp).fill(
			    pro.p(), (pred_pro_p_ang - pro.p()))
			    //okay
			    histos.computeIfAbsent('p_pro_dpp_pro_from_beam_eangle_' + ctof + '_' + sector_pro + '_' + isr_fsr_flag, histoBuilders2.p_pro_dp).fill(
			    pro.p(), (pred_pro_p_ang - pro.p())/pro.p())

			    //fill 1d histograms for fitting later
			    //okay
			    histos.computeIfAbsent('dp_pro_from_angles_' + ctof + '_' + sector_pro + '_thetabin' + theta_bin_pro + '_' + isr_fsr_flag, histoBuilders.dp_pro).fill((pred_pro_p - pro.p())/pro.p())
			    //okay
			    //histos.computeIfAbsent('dp_pro_from_angles_' + ctof + '_' + sector_pro + '_thetabin' + theta_bin_pro + '_phibin' + phi_bin_pro + '_' + isr_fsr_flag, histoBuilders.dp_pro).fill((pred_pro_p - pro.p())/pro.p())

			    // David - only use the _from_angles_ histograms
			    //prob not okay
 	 		    histos.computeIfAbsent('p_pro_dp_pro_from_proangles_' + ctof + '_' + sector_pro + '_' + isr_fsr_flag, histoBuilders2.p_pro_dp).fill(
			    pro.p(), (pred_pro_p_beam_prang - pro.p()))
 			    histos.computeIfAbsent('p_pro_dpp_pro_from_proangles_' + ctof + '_' + sector_pro + '_' + isr_fsr_flag, histoBuilders2.p_pro_dp).fill(
			    pro.p(), (pred_pro_p_beam_prang - pro.p())/pro.p() )

			    // Check deltap/p vs p in different theta bins of proton
 			    histos.computeIfAbsent('p_pro_dpp_pro_from_angles_' + ctof + '_' + sector_pro + '_thetabin' + theta_bin_pro + '_' + isr_fsr_flag, histoBuilders2.p_pro_dp).fill(
			    pro.p(), (pred_pro_p_beam_prang - pro.p())/pro.p() )

			    //not okay bc uses electron momentum which is altered via radiation in final state, and not correct with ISR
			    //histos.computeIfAbsent('theta_pro_dtheta_pro_' + ctof + '_' + sector + '_' + isr_fsr_flag, histoBuilders2.theta_pro_dtheta).fill(
			    //	Math.toDegrees(pro.theta()), Math.toDegrees(pred_pro.theta - pro.theta())/Math.toDegrees(pro.theta()))

			   //okay
			   histos.computeIfAbsent('theta_pro_dtheta_pro_from_beam_eangle_' + ctof + '_' + sector_pro + '_' + isr_fsr_flag, histoBuilders2.theta_pro_dtheta).fill(
				Math.toDegrees(pro.theta()), Math.toDegrees(pred_pro_theta_ang - pro.theta()))
			   //okay
			   histos.computeIfAbsent('theta_pro_dthetatheta_pro_from_beam_eangle_' + ctof + '_' + sector_pro + '_' + isr_fsr_flag, histoBuilders2.theta_pro_dtheta).fill(
				Math.toDegrees(pro.theta()), Math.toDegrees(pred_pro_theta_ang - pro.theta())/Math.toDegrees(pro.theta()))
			    
			   histos.computeIfAbsent('chi2_dp_ele_' + ctof + '_' + isr_fsr_flag, histoBuilders2.chi2_dp_ele).fill(event.chi2pid[it], pred_ele_p - ele.p())
			   histos.computeIfAbsent('chi2_dp_pro_' + ctof + '_' + isr_fsr_flag, histoBuilders2.chi2_dp_ele).fill(event.chi2pid[it], pred_pro_p - pro.p())

			    

			}
		    }
		    		   

		    //////////////////////////////////////////////////////////////////////////
		    //////////////////////////////////////////////////////////////////////////
 		    // Everything Else Passed (eep) 		    
		    //These are radiated elastic events - before selecting the ISR or FSR
		    if( isr_fsr_flag != null && pass_angle_ep && (pass_missing_mass_rad || pass_missing_mass ) ){
			histos.computeIfAbsent('w_eep_rad_elastic_' + ctof, histoBuilders.w).fill(pkin.w)
		    }
		    if( isr_fsr_flag != null && pass_angle_ep && pass_high_w){
			// define missing mass from pkin for fsr and effective_missing_mass for isr
			def corr_mm2 = pkin.missing_mass
			if( isr_fsr_flag == 'isr' ) corr_mm2 = effective_pkin.missing_mass
			histos.computeIfAbsent('missing_mass_eep_rad_elastic_' + ctof, histoBuilders.missing_mass_small).fill(corr_mm2)//pkin.missing_mass)
		    }
		    if( isr_fsr_flag != null && pass_missing_mass_rad && pass_high_w ){
			histos.computeIfAbsent('angle_ep_eep_rad_elastic_' + ctof, histoBuilders.angle_ep).fill(pkin.angle)
		    }
 		    if( isr_fsr_flag != null && pass_missing_mass_rad ){
			histos.computeIfAbsent('angle_ep_eep_rad_elastic_nowcut_' + ctof, histoBuilders.angle_ep).fill(pkin.angle)
		    }
		    if( isr_fsr_flag != null && pass_angle_ep ){
			histos.computeIfAbsent('angle_ep_eep_rad_elastic_noangleepcut__' + ctof, histoBuilders.angle_ep).fill(pkin.angle)
		    }

		    if( pass_angle_ep ){
			def corr_mm2=0
			if( isr_fsr_flag == 'isr') corr_mm2 = effective_pkin.missing_mass
			else if( isr_fsr_flag == 'fsr' ) corr_mm2 = pkin.missing_mass
			else corr_mm2 = pkin.missing_mass
			histos.computeIfAbsent('missing_mass_eep_rad_elastic_pass_angleep_' + ctof, histoBuilders.missing_mass_small).fill(cor_mm2)//pkin.missing_mass)
			histos.computeIfAbsent('w_rad_elastic_pass_angle_ep_' + ctof, histoBuilders.w).fill(pkin.w)
		    }
		    if ( pass_high_w ){
			histos.computeIfAbsent('angle_ep_rad_elastic_pass_high_w_' + ctof, histoBuilders.angle_ep).fill(pkin.angle)
		    }
		    if( pass_angle_ep && pass_high_w ){
			def corr_mm2=0
			if( isr_fsr_flag == 'isr') corr_mm2 = effective_pkin.missing_mass
			else if( isr_fsr_flag == 'fsr' ) corr_mm2 = pkin.missing_mass
			else corr_mm2 = pkin.missing_mass
 			histos.computeIfAbsent('missing_mass_eep_rad_elastic_pass_angleep_highw_' + ctof, histoBuilders.missing_mass_small).fill(cor_mm2)//pkin.missing_mass)
			histos.computeIfAbsent('missing_mass_eep_rad_elastic_pass_angleep_highw_bad' + ctof, histoBuilders.missing_mass_small).fill(pkin.missing_mass)
		    }

		    if( pass_missing_mass_rad || pass_missing_mass ){
			histos.computeIfAbsent('angle_ep_eep_' + ctof + '_pass_both_missing_mass', histoBuilders.angle_ep).fill(pkin.angle)
		    }
			
		    		    
		    //ISR
		    if (pass_angle_ep && pass_missing_mass_rad && pass_theta_gamma){ //remove
			histos.computeIfAbsent('w_eep_' + ctof, histoBuilders.w).fill(pkin.w)
		    }

		    if (pass_missing_mass_rad && pass_high_w && pass_theta_gamma){
			histos.computeIfAbsent('angle_ep_eep_' + ctof, histoBuilders.angle_ep).fill(pkin.angle)
		    }

		    if ( pass_angle_ep && pass_high_w ){
			histos.computeIfAbsent('theta_gamma_eep_' + ctof, histoBuilders.theta_gamma).fill(pkin.theta_gamma)
			if (pass_missing_mass_rad || pass_missing_mass ){
 			    histos.computeIfAbsent('w_theta_sum_eepb_theta_gamma' + ctof, histoBuilders2.w_theta_sum).fill(
				pkin.w, Math.toDegrees(ele.theta() + pro.theta())
			    }
			)
		    }
		    if (pass_high_w && pass_theta_gamma && pass_angle_ep){
			histos.computeIfAbsent('missing_mass_eep_' + ctof, histoBuilders.missing_mass_small).fill(cor_mm2)
		    }

		    //FSR
		    if (pass_angle_ep && pass_missing_mass_rad && pass_theta_egamma){
			histos.computeIfAbsent('w_eep_' + ctof + '_fsr', histoBuilders.w).fill(pkin.w)
		    }

 		    if (pass_missing_mass_rad && pass_high_w && pass_theta_egamma){
			histos.computeIfAbsent('angle_ep_eep_' + ctof  + '_fsr', histoBuilders.angle_ep).fill(pkin.angle)
		    }

		    if ( (pass_missing_mass_rad || pass_missing_mass ) && pass_angle_ep && pass_high_w){
			histos.computeIfAbsent('theta_egamma_eep_' + ctof, histoBuilders.theta_gamma).fill(pkin.theta_egamma)			
		    }		    

		    if (pass_high_w && pass_theta_egamma && pass_angle_ep){						histos.computeIfAbsent('missing_mass_eep_' + ctof  + '_fsr', histoBuilders.missing_mass_small).fill(cor_mm2)
		    }

		    /////////////////////////////////////////
		    // Everything else pass but for elastic
		    if (pass_angle_ep && pass_missing_mass ){
			histos.computeIfAbsent('w_eep_' + ctof + '_elastic', histoBuilders.w).fill(pkin.w)
			histos.computeIfAbsent('theta_egamma_eep_' + ctof + '_pass_aep_mm2', histoBuilders.theta_gamma).fill(pkin.theta_egamma)			
			histos.computeIfAbsent('theta_gamma_eep_' + ctof + '_pass_aep_mm2', histoBuilders.theta_gamma).fill(pkin.theta_gamma)			
		    }

 		    if (pass_missing_mass && pass_w_elastic){
			histos.computeIfAbsent('angle_ep_eep_' + ctof  + '_elastic', histoBuilders.angle_ep).fill(pkin.angle)
		    }

 		    if (pass_missing_mass && pass_angle_ep && pass_w_elastic){
			histos.computeIfAbsent('theta_egamma_eep_' + ctof + '_elastic', histoBuilders.theta_gamma).fill(pkin.theta_egamma)			
			histos.computeIfAbsent('theta_gamma_eep_' + ctof + '_elastic', histoBuilders.theta_gamma).fill(pkin.theta_gamma)			
		    }		    
		    if (pass_w_elastic && pass_angle_ep){
			histos.computeIfAbsent('missing_mass_eep_' + ctof  + '_elastic', histoBuilders.missing_mass_small).fill(pkin.missing_mass)
		    }
		    

		    
		    
                }
            }
            eventIndex++
        }
    }
}

def out = new TDirectory()
out.mkdir("histos")
out.cd("histos")
histos.values().each { out.addDataSet(it) }
out.writeFile("mon-kcor.hipo")

def root_out = new ROOTFile("mon-kcor-rga18v5.root")
histos.each{root_out.writeDataSet(it.value)}
root_out.close()
