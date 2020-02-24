package mon
import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import java.util.concurrent.ConcurrentHashMap
import java.util.LinkedHashMap
import pid.electron.ElectronFromEvent
import event.Event
import event.EventConverter
import utils.KinTool

class ElectronMon{

    def beam = 10.604

    //call electron cut constructor
    def electron = new ElectronFromEvent();

    def myElectronCutStrategies = [
	electron.passElectronStatus,
	electron.passElectronChargeCut,
	electron.passElectronEBPIDCut,
	electron.passElectronNpheCut,
	electron.passElectronVertexCut,
	electron.passElectronPCALFiducialCut,
	electron.passElectronEIEOCut,
	electron.passElectronDCR1,
	electron.passElectronDCR2,
	electron.passElectronDCR3
    ]

    // Combinations of Cut Pass Configurations to Plot
    // charge, status and pass DCr1, r2, r3 fid cuts
    def pass_dc_comb=[true,true,true,true,false,false,false,true,true,true]
    // charge, status and pass PCAL fid DCr1, r2, r3 fid cuts
    def pass_dc_pcal_comb=[true,true,true,true,false,true,false,true,true,true]
    def pass_cut_comb = [0:pass_dc_comb,
			 1:pass_dc_pcal_comb]
    
    def histos = new ConcurrentHashMap()
    def histoBuilders = [

	ptheta : {title -> new H2F("$title","$title", 900, 0, beam, 900, 0.0, 65.0 ) },
	phitheta : {title -> new H2F("$title","$title", 900, 180, -180, 900, 0.0, 65.0 ) },
	vzphi : { title -> new H2F("$title","$title", 600, -180, 180, 600, -20, 50 ) },
	vztheta : { title -> new H2F("$title","$title", 1000, -100, 100, 900, 0.0, 65.0 ) },
	hitxy : { title -> new H2F("$title","$title", 900, -450.0, 450.0, 900, -450.0, 450.0 ) },
	ecsf : { title -> new H2F("$title","$title", 500, 0.0, beam, 500, 0.0, 0.3) },
	uvwsf : { title -> new H2F("$title","$title", 500, 0.0, 500.0, 500, 0.0, 0.3) },
	tchi2 : { title -> new H1F("$title","$title", 500, 0.0, 500 ) }, 
	tndf : { title -> new H1F("$title","$title", 50, 0.0, 50 ) },
	nphe : { title -> new H1F("$title","$title", 40, 0.0, 40.0 ) },
	eieo : { title -> new H2F("$title","$title", 1000, 0.0, 1.0, 1000, 0.0, 1.0 ) },
	vz   : { title -> new H1F("$title","$title", 400, -40, 40) },
	vxy : { title -> new H2F("$title","$title", 400, -1, 1,400, -1, 1) },

	time : { title -> new H1F("$title","$title", 2000, -450, 450) },
	path : { title -> new H1F("$title","$title", 2000, -400, 400) },//
	sectors : { title -> new H1F("$title","$title", 7, 0, 7) },
	helicity : { title -> new H1F("$title","$title", 6, -3 , 3) },
	w    : { title -> new H1F("$title","$title", 500, 0.5, beam+0.3 ) },
	q2w  : { title -> new H2F("$title","$title", 500, 0.5, beam+0.3, 500, 0, 2 * beam ) },  
	passrate : {title -> new H1F("$title","$title", 8, 0.0, 8.0 ) },
	wtheta : { title -> new H2F("$title","$title", 100, 0.0, 2.5, 100, 0.0, 20.0 ) },
	deltakin : {title -> new H1F("$title","$title", 500, -2, 2  ) },
	deltakin2D : {title -> new H2F("$title","$title", 1000, 0.0, 500, 500, -2, 2  ) },
	deltakin2D_beam : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, -2, 2  ) }

    ]

    def processEvent( event ){
	
	def my_el_cuts = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()	  
	
	my_el_cuts.each{ index, value ->
	    def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
	    def p = lv.mag()
	    def vz = event.vz[index]
	    def theta = Math.toDegrees(lv.theta())
	    def phi = Math.toDegrees(lv.phi())
	    
	    value.eachWithIndex{ cut_res, cut_ind -> cut_res ? histos.computeIfAbsent("passResults",histoBuilders.passrate).fill(cut_ind) : null }

	    if (event.cherenkov_status.contains(index)) { 
		def nphe = event.nphe[index]	    
		value.eachWithIndex{ cut_result, cut_index ->
		    if(cut_result){histos.computeIfAbsent("nphe_cut$cut_index",histoBuilders.nphe).fill(nphe)}
		}
		//plot nphe for tracks passing all cuts
		if(!value.contains(false)){histos.computeIfAbsent("nphe_passall",histoBuilders.nphe).fill(nphe)}		    
	    }
	    
	    
	    if( event.pcal_status.contains(index) ){
		def e_pcal = event.pcal_energy[index]
		def pcal_sect = event.pcal_sector[index]
		def pcal_u = event.pcal_u[index]
		def pcal_v = event.pcal_v[index]
		def pcal_w = event.pcal_w[index]
		def pcal_t = event.pcal_time[index]
		def pcal_l = event.pcal_path[index]
		value.eachWithIndex{ cut_result, cut_index ->
		    if(cut_result){
			histos.computeIfAbsent("pcalsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.ecsf).fill(p, e_pcal/p) 
			histos.computeIfAbsent("pcalusf_cut$cut_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_u, e_pcal/p)
			histos.computeIfAbsent("pcalvsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_v, e_pcal/p)
			histos.computeIfAbsent("pcalwsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_w, e_pcal/p)
			histos.computeIfAbsent("pcalut_cut$cut_index"+"_s$pcal_sect",histoBuilders.time).fill(pcal_u, pcal_t)
			histos.computeIfAbsent("pcalvt_cut$cut_index"+"_s$pcal_sect",histoBuilders.time).fill(pcal_v, pcal_t)
			histos.computeIfAbsent("pcalwt_cut$cut_index"+"_s$pcal_sect",histoBuilders.time).fill(pcal_w, pcal_t)
		    }
		}
		if(!value.contains(false)){
		    histos.computeIfAbsent("pcalsf_passall"+"_s$pcal_sect",histoBuilders.ecsf).fill(p, e_pcal/p)
		    histos.computeIfAbsent("pcalusf_passall"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_u, e_pcal/p)
		    histos.computeIfAbsent("pcalvsf_passall"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_v, e_pcal/p)
		    histos.computeIfAbsent("pcalwsf_passall"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_w, e_pcal/p)

		    histos.computeIfAbsent("pcalut_passall"+"_s$pcal_sect",histoBuilders.time).fill(pcal_u, pcal_t)
		    histos.computeIfAbsent("pcalvt_passall"+"_s$pcal_sect",histoBuilders.time).fill(pcal_v, pcal_t)
		    histos.computeIfAbsent("pcalwt_passall"+"_s$pcal_sect",histoBuilders.time).fill(pcal_w, pcal_t)
		}		
		pass_cut_comb.each{ cut_comb_index, cut_comb_value ->
		    if( cut_comb_value.equals(value) ){
			histos.computeIfAbsent("pcalsf_cut_comb$cut_comb_index"+"_s$pcal_sect",histoBuilders.ecsf).fill(p, e_pcal/p) 
			histos.computeIfAbsent("pcalusf_cut_comb$cut_comb_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_u, e_pcal/p)
			histos.computeIfAbsent("pcalvsf_cut_comb$cut_comb_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_v, e_pcal/p)
			histos.computeIfAbsent("pcalwsf_cut_comb$cut_comb_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_w, e_pcal/p)
			histos.computeIfAbsent("pcalut_cut_comb$cut_comb_index"+"_s$pcal_sect",histoBuilders.time).fill(pcal_u, pcal_t)
			histos.computeIfAbsent("pcalvt_cut_comb$cut_comb_index"+"_s$pcal_sect",histoBuilders.time).fill(pcal_v, pcal_t)
			histos.computeIfAbsent("pcalwt_cut_comb$cut_comb_index"+"_s$pcal_sect",histoBuilders.time).fill(pcal_w, pcal_t)
		    }
		}

	    }

	    if( event.ecal_inner_status.contains(index) &&  event.ecal_outer_status.contains(index)){
		def e_eical = event.ecal_inner_energy[index]
 		def e_eocal = event.ecal_outer_energy[index]
		def ei_sect = event.ecal_inner_sector[index]
		def eo_sect = event.ecal_outer_sector[index]
		
		if( ei_sect == eo_sect ){			
		    value.eachWithIndex{ cut_result, cut_index ->
			if(cut_result){histos.computeIfAbsent("eieo_cut$cut_index"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)}
		    }
		    if(!value.contains(false)){histos.computeIfAbsent("eieo_passall"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)}		    
		}		    
	    }
	    
	    if( event.ecal_inner_status.contains(index) ){
		def e_ecal = event.ecal_inner_energy[index]
		def eical_sect = event.ecal_inner_sector[index]
		def eical_u = event.ecal_inner_u[index]
		def eical_v = event.ecal_inner_v[index]
		def eical_w = event.ecal_inner_w[index]
		def eical_t = event.ecal_inner_time[index]
		def eical_l = event.ecal_inner_path[index]
		value.eachWithIndex{ cut_result, cut_index ->
		    if(cut_result){
			histos.computeIfAbsent("eicalsf_cut$cut_index"+"_s$eical_sect",histoBuilders.ecsf).fill(p, e_ecal/p) 
			histos.computeIfAbsent("eical_usf_cut$cut_index"+"_s$eical_sect",histoBuilders.uvwsf).fill(eical_u, e_ecal/p)
			histos.computeIfAbsent("eical_vsf_cut$cut_index"+"_s$eical_sect",histoBuilders.uvwsf).fill(eical_v, e_ecal/p)
			histos.computeIfAbsent("eical_wsf_cut$cut_index"+"_s$eical_sect",histoBuilders.uvwsf).fill(eical_w, e_ecal/p)

			histos.computeIfAbsent("eical_ut_cut$cut_index"+"_s$eical_sect",histoBuilders.time).fill(eical_u, eical_t)
			histos.computeIfAbsent("eical_vt_cut$cut_index"+"_s$eical_sect",histoBuilders.time).fill(eical_v, eical_t)
			histos.computeIfAbsent("eical_wt_cut$cut_index"+"_s$eical_sect",histoBuilders.time).fill(eical_w, eical_t)
		    }			    
		}
		if(!value.contains(false)){
		    histos.computeIfAbsent("eicalsf_passall_s$eical_sect",histoBuilders.ecsf).fill(p, e_ecal/p) 
		    histos.computeIfAbsent("eical_usf_passall_s$eical_sect",histoBuilders.uvwsf).fill(eical_u, e_ecal/p)
		    histos.computeIfAbsent("eical_vsf_passall_s$eical_sect",histoBuilders.uvwsf).fill(eical_v, e_ecal/p)
		    histos.computeIfAbsent("eical_wsf_passall_s$eical_sect",histoBuilders.uvwsf).fill(eical_w, e_ecal/p)
		    
		    histos.computeIfAbsent("eical_ut_passall_s$eical_sect",histoBuilders.time).fill(eical_u, eical_t)
		    histos.computeIfAbsent("eical_vt_passall_s$eical_sect",histoBuilders.time).fill(eical_v, eical_t)
		    histos.computeIfAbsent("eical_wt_passall_s$eical_sect",histoBuilders.time).fill(eical_w, eical_t)
		}

		pass_cut_comb.each{ cut_comb_index, cut_comb_value ->
		    if( cut_comb_value.equals(value) ){
			histos.computeIfAbsent("eicalsf_cut_comb$cut_comb_index"+"_s$eical_sect",histoBuilders.ecsf).fill(p, e_ecal/p) 
			histos.computeIfAbsent("eical_usf_cut_comb$cut_comb_index"+"_s$eical_sect",histoBuilders.uvwsf).fill(eical_u, e_ecal/p)
			histos.computeIfAbsent("eical_vsf_cut_comb$cut_comb_index"+"_s$eical_sect",histoBuilders.uvwsf).fill(eical_v, e_ecal/p)
			histos.computeIfAbsent("eical_wsf_cut_comb$cut_comb_index"+"_s$eical_sect",histoBuilders.uvwsf).fill(eical_w, e_ecal/p)
			
			histos.computeIfAbsent("eical_ut_cut_comb$cut_comb_index"+"_s$eical_sect",histoBuilders.time).fill(eical_u, eical_t)
			histos.computeIfAbsent("eical_vt_cut_comb$cut_comb_index"+"_s$eical_sect",histoBuilders.time).fill(eical_v, eical_t)
			histos.computeIfAbsent("eical_wt_cut_comb$cut_comb_index"+"_s$eical_sect",histoBuilders.time).fill(eical_w, eical_t)
		    }
		}
	    }
	    
	    if( event.ecal_outer_status.contains(index) ){
		def e_ecal = event.ecal_outer_energy[index]
		def eocal_sect = event.ecal_outer_sector[index]
		def eocal_u = event.ecal_outer_u[index]
		def eocal_v = event.ecal_outer_v[index]
		def eocal_w = event.ecal_outer_w[index]
		def eocal_t = event.ecal_outer_time[index]
		def eocal_l = event.ecal_outer_path[index]
		value.eachWithIndex{ cut_result, cut_index ->
		    if(cut_result){
			histos.computeIfAbsent("eocalsf_cut$cut_index"+"_s$eocal_sect",histoBuilders.ecsf).fill(p, e_ecal/p) 
			histos.computeIfAbsent("eocal_usf_cut$cut_index"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_u, e_ecal/p)
			histos.computeIfAbsent("eocal_vsf_cut$cut_index"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_v, e_ecal/p)
			histos.computeIfAbsent("eocal_wsf_cut$cut_index"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_w, e_ecal/p)
			
			histos.computeIfAbsent("eocal_ut_cut$cut_index"+"_s$eocal_sect",histoBuilders.time).fill(eocal_u, eocal_t)
			histos.computeIfAbsent("eocal_vt_cut$cut_index"+"_s$eocal_sect",histoBuilders.time).fill(eocal_v, eocal_t)
			histos.computeIfAbsent("eocal_wt_cut$cut_index"+"_s$eocal_sect",histoBuilders.time).fill(eocal_w, eocal_t)

		    }
		}
		if(!value.contains(false)){
		    histos.computeIfAbsent("eocalsf_passall_s$eocal_sect",histoBuilders.ecsf).fill(p, e_ecal/p) 
		    histos.computeIfAbsent("eocal_usf_passall"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_u, e_ecal/p)
		    histos.computeIfAbsent("eocal_vsf_passall"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_v, e_ecal/p)
		    histos.computeIfAbsent("eocal_wsf_passall"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_w, e_ecal/p)
		    
		    histos.computeIfAbsent("eocal_ut_passall"+"_s$eocal_sect",histoBuilders.time).fill(eocal_u, eocal_t)
		    histos.computeIfAbsent("eocal_vt_passall"+"_s$eocal_sect",histoBuilders.time).fill(eocal_v, eocal_t)
		    histos.computeIfAbsent("eocal_wt_passall"+"_s$eocal_sect",histoBuilders.time).fill(eocal_w, eocal_t)
		}

		pass_cut_comb.each{ cut_comb_index, cut_comb_value ->
		    if( cut_comb_value.equals(value) ){
			histos.computeIfAbsent("eocalsf_cut_comb$cut_comb_index"+"_s$eocal_sect",histoBuilders.ecsf).fill(p, e_ecal/p) 
			histos.computeIfAbsent("eocal_usf_cut_comb$cut_comb_index"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_u, e_ecal/p)
			histos.computeIfAbsent("eocal_vsf_cut_comb$cut_comb_index"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_v, e_ecal/p)
			histos.computeIfAbsent("eocal_wsf_cut_comb$cut_comb_index"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_w, e_ecal/p)

			histos.computeIfAbsent("eocal_ut_cut_comb$cut_comb_index"+"_s$eocal_sect",histoBuilders.time).fill(eocal_u, eocal_t)
			histos.computeIfAbsent("eocal_vt_cut_comb$cut_comb_index"+"_s$eocal_sect",histoBuilders.time).fill(eocal_v, eocal_t)
			histos.computeIfAbsent("eocal_wt_cut_comb$cut_comb_index"+"_s$eocal_sect",histoBuilders.time).fill(eocal_w, eocal_t)
		    }
		}
	    }

	    if( event.dc1_status.contains(index) ){		    
		def hits = event.dc1[index]
		def hit = hits.find{it.layer==12}.each{			
		    value.eachWithIndex{ cut_result, cut_index ->
			if(cut_result){
 			    if( cut_result ){ histos.computeIfAbsent("dcr1hit_cut$cut_index",histoBuilders.hitxy).fill(it.x, it.y) }
			}
		    }		    
		    if(!value.contains(false)){histos.computeIfAbsent("dcr1hit_passall",histoBuilders.hitxy).fill(it.x, it.y)}
		    
		    
		    pass_cut_comb.each{ cut_comb_index, cut_comb_value ->
			if( cut_comb_value.equals(value) ){
			    histos.computeIfAbsent("dcr1hit_cut_comb$cut_comb_index",histoBuilders.hitxy).fill(it.x, it.y)
			}
		    }
		}

	    }
 	    
	    if( event.dc2_status.contains(index) ){		    
		def hits = event.dc1[index]
		def hit = hits.find{it.layer==24}.each{
		    value.eachWithIndex{ cut_result, cut_index ->
			if(cut_result){
			    if( cut_result ){ histos.computeIfAbsent("dcr2hit_cut$cut_index",histoBuilders.hitxy).fill(it.x, it.y)}
			}
			if(!value.contains(false)){histos.computeIfAbsent("dcr2hit_passall",histoBuilders.hitxy).fill(it.x, it.y)}
		    }
		
		    pass_cut_comb.each{ cut_comb_index, cut_comb_value ->
			if( cut_comb_value.equals(value) ){
			    histos.computeIfAbsent("dcr2hit_cut_comb$cut_comb_index",histoBuilders.hitxy).fill(it.x, it.y)
			}
		    }
		}
	    }

	    if( event.dc3_status.contains(index) ){		    
		def hits = event.dc1[index]
		def hit = hits.find{it.layer==36}.each{
		    value.eachWithIndex{ cut_result, cut_index ->
			if(cut_result){
			    if( cut_result ){ histos.computeIfAbsent("dcr3hit_cut$cut_index",histoBuilders.hitxy).fill(it.x, it.y) }
			}
			if(!value.contains(false)){histos.computeIfAbsent("dcr3hit_passall",histoBuilders.hitxy).fill(it.x, it.y)}
		    }
		    pass_cut_comb.each{ cut_comb_index, cut_comb_value ->
			if( cut_comb_value.equals(value) ){
			    histos.computeIfAbsent("dcr3hit_cut_comb$cut_comb_index",histoBuilders.hitxy).fill(it.x, it.y)
			}
		    }
		}
	    }


	    if( event.tof_status.contains(index) ){
		def ftof_sector = event.tof_sector[index]
		def ftof_t = event.tof_time[index]
		def ftof_p = event.tof_path[index]
		def ftof_layer = event.tof_paddle[index]
		def ftof_t_corr = ftof_t - ftof_p/(29.9792458)
		value.eachWithIndex{ cut_result, cut_index ->
		    if(cut_result){		   
			histos.computeIfAbsent("ftoft_cut$cut_index"+"_s$ftof_sector",histoBuilders.time).fill(ftof_t)
			histos.computeIfAbsent("ftoft_cut$cut_index"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.time).fill(ftof_t)
			histos.computeIfAbsent("ftof_t_corr_cut$cut_index"+"_s$ftof_sector",histoBuilders.time).fill(ftof_t_corr)
 			histos.computeIfAbsent("ftof_t_corr_cut$cut_index"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.time).fill(ftof_t_corr)
		    }
		}
		if(!value.contains(false)){
		    histos.computeIfAbsent("ftoft_passall_s$ftof_sector",histoBuilders.time).fill(ftof_t)
		    histos.computeIfAbsent("ftoft_passall_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.time).fill(ftof_t)		    
		    histos.computeIfAbsent("ftof_t_corr_passall_s$ftof_sector",histoBuilders.time).fill(ftof_t_corr)
 		    histos.computeIfAbsent("ftof_t_corr_passall_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.time).fill(ftof_t_corr)
		}
		pass_cut_comb.each{ cut_comb_index, cut_comb_value ->
 		    if( cut_comb_value.equals(value) ){
			histos.computeIfAbsent("ftoft_cut_comb$cut_comb_index"+"_s$ftof_sector",histoBuilders.time).fill(ftof_t)
			histos.computeIfAbsent("ftoft_cut_comb$cut_comb_index"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.time).fill(ftof_t)
			histos.computeIfAbsent("ftof_t_corr_cut_comb$cut_comb_index"+"_s$ftof_sector",histoBuilders.time).fill(ftof_t_corr)
 			histos.computeIfAbsent("ftof_t_corr_cut_comb$cut_comb_index"+"_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.time).fill(ftof_t_corr)
		    }
		}
	    }
	    
	    	    	    
	    //generic p, theta, phi, vx, vy, vz information
	    //also CHI2 and NDF
	    if( event.npart>0 && event.dc1_status.contains(index)) {
		def dc_sector  = event.dc_sector[index]		    
		def dc_chi2 = event.dc_chi2[index]
		def dc_ndf = event.dc_ndf[index]
		def vx = event.vx[index]
		def vy = event.vy[index]
		def measured_el = LorentzVector.withPID(11,event.px[index],event.py[index],event.pz[index])
		def delta_energy = KinTool.delta_meas_energy(beam, measured_el) 
		def delta_theta = KinTool.delta_meas_theta(beam, measured_el)
		
		value.eachWithIndex{ cut_result, cut_index ->
		    histos.computeIfAbsent("vz_cut$cut_index"+"_s$dc_sector",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("ptheta_cut$cut_index"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
		    histos.computeIfAbsent("phitheta_cut$cut_index",histoBuilders.phitheta).fill(phi, theta)
		    histos.computeIfAbsent("vztheta_cut$cut_index"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vxy_cut$cut_index"+"_s$dc_sector",histoBuilders.vxy).fill(vx,vy)
		    histos.computeIfAbsent("vzphi_cut$cut_index"+"_s$dc_sector",histoBuilders.vzphi).fill(phi,vz)
		    histos.computeIfAbsent("chi_cut$cut_index"+"_s$dc_sector",histoBuilders.tchi2).fill(dc_chi2)
 		    histos.computeIfAbsent("ndf_cut$cut_index"+"_s$dc_sector",histoBuilders.tndf).fill(dc_ndf)		    
		}	    	

		if(!value.contains(false)){		    
		    histos.computeIfAbsent("vz_passall"+"_s$dc_sector",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("ptheta_passall"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
		    histos.computeIfAbsent("phitheta_passall",histoBuilders.phitheta).fill(phi, theta)
		    histos.computeIfAbsent("vztheta_passall"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vzphi_passall",histoBuilders.vzphi).fill(phi,vz)
 		    histos.computeIfAbsent("chi_passall"+"_s$dc_sector",histoBuilders.tchi2).fill(dc_chi2)
 		    histos.computeIfAbsent("ndf_passall"+"_s$dc_sector",histoBuilders.tndf).fill(dc_ndf)		
 		    histos.computeIfAbsent("delta_energy_passall_$dc_sector",histoBuilders.deltakin).fill(delta_energy)
		    histos.computeIfAbsent("delta_theta_passall_$dc_sector",histoBuilders.deltakin).fill(delta_theta)
 		    histos.computeIfAbsent("delta_energy_vs_chi2_passall_$dc_sector",histoBuilders.deltakin2D).fill(dc_chi2,delta_energy)
		    histos.computeIfAbsent("delta_theta_vs_chi2_passall_$dc_sector",histoBuilders.deltakin2D).fill(dc_chi2,delta_theta)
		    histos.computeIfAbsent("delta_energy_vs_p_passall_$dc_sector",histoBuilders.deltakin2D_beam).fill(p,delta_energy)
		    histos.computeIfAbsent("delta_theta_vs_p_passall_$dc_sector",histoBuilders.deltakin2D_beam).fill(p,delta_theta)
		}

		pass_cut_comb.each{ cut_comb_index, cut_comb_value ->
 		    if( cut_comb_value.equals(value) ){
			histos.computeIfAbsent("vz_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.vz).fill(vz)
			histos.computeIfAbsent("ptheta_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
			histos.computeIfAbsent("phitheta_cut_comb$cut_comb_index",histoBuilders.phitheta).fill(phi, theta)
			histos.computeIfAbsent("vztheta_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
			histos.computeIfAbsent("vxy_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.vxy).fill(vx,vy)
			histos.computeIfAbsent("vzphi_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.vzphi).fill(phi,vz)
			histos.computeIfAbsent("chi_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.tchi2).fill(dc_chi2)
 			histos.computeIfAbsent("ndf_cut_comb$cut_comb_index"+"_s$dc_sector",histoBuilders.tndf).fill(dc_ndf)		    
			histos.computeIfAbsent("delta_energy_cut_comb$cut_comb_index"+"_$dc_sector",histoBuilders.deltakin).fill(delta_energy)
			histos.computeIfAbsent("delta_theta_cut_comb$cut_comb_index"+"_$dc_sector",histoBuilders.deltakin).fill(delta_theta)
 			histos.computeIfAbsent("delta_energy_vs_chi2_cut_comb$cut_comb_index"+"_$dc_sector",histoBuilders.deltakin2D).fill(dc_chi2,delta_energy)
			histos.computeIfAbsent("delta_theta_vs_chi2_cut_comb$cut_comb_index"+"_$dc_sector",histoBuilders.deltakin2D).fill(dc_chi2,delta_theta)
			histos.computeIfAbsent("delta_energy_vs_p_cut_comb$cut_comb_index"+"_$dc_sector",histoBuilders.deltakin2D_beam).fill(p,delta_energy)
			histos.computeIfAbsent("delta_theta_vs_p_cut_comb$cut_comb_index"+"_$dc_sector",histoBuilders.deltakin2D_beam).fill(p,delta_theta)
		    }
		}
	    }	
	}

	if(event.npart>0){		
	    histos.computeIfAbsent("event_helicity",histoBuilders.helicity).fill(event.helicity)
	    histos.computeIfAbsent("event_start_time",histoBuilders.time).fill(event.start_time)
	    histos.computeIfAbsent("event_rf_time",histoBuilders.time).fill(event.rf_time)		    		
	    event.charge.findResults{ it.value < 0 ? it.value : null}?.each{histos.computeIfAbsent("event_charge_neg",histoBuilders.helicity).fill(it)}
	    event.charge.findResults{ it.value > 0 ? it.value : null}?.each{histos.computeIfAbsent("event_charge_pos",histoBuilders.helicity).fill(it)}
	}		
    }

}
