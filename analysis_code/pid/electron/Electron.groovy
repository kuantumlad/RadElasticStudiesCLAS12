package pid.electron
import org.jlab.detector.base.DetectorType

class Electron {
    
    def ebPID = 11
    def min_nphe = 2
    
    //vertex wide enough for all sectors
    def min_vz = -12
    def max_vz = 9
    
    def min_pcal_dep = 0.06

    //calorimeter fiducial cut
    def min_u = 0
    def max_u = 420
    def min_v = 0
    def max_v = 420
    def min_w = 0
    def max_w = 420

    //dcr1,2,3 fiducial cut
    //require functional form
    def sect_angle_coverage = 60
    def heightR1 = 15 // height is lower for outbending runs - this parameter is field config. / sector dependent
    def radiusR1 = 25 

    def heightR2 = 35
    def radiusR2 = 40

    def heightR3 = 48
    def radiusR3 = 49
    
    static def findElectron = { event ->
	def pbank = event.getBank("REC::Particle")
	return (0..<pbank.rows())
	    .find{pbank.getInt('pid',it)==11 && pbank.getShort('status',it)<0}
    }

    def passElectronStatus = { bank, index -> 
	return (bank.part.getInt('status',index)<0 )
    }

    def passElectronEBPIDCut = {bank, index ->
	return (bank.part.getInt('pid',index) == ebPID)
    }

    def passElectronChargeCut= {bank, index ->
	return (bank.part.getByte('charge',index) < 0)
    }

    def passElectronNpheCut= {bank, index ->
	return (0..<bank.cc.rows()).any{bank.cc.getInt('pindex',it) == index &&
					bank.cc.getFloat('nphe',it) > min_nphe}
    }

    def passElectronVertexCut= {bank, index ->
	return (bank.part.getFloat('vz',index).with{it < max_vz && it > min_vz})
    }

    def passElectronPCALFiducialCut= {bank, index ->		
	// can probably be changed to a find like closure to avoid looping over all indices
	return (0..<bank.ec.rows()).any{bank.ec.getByte('detector',it) == DetectorType.ECAL.getDetectorId() &&
					bank.ec.getInt('layer',it) == 1 &&
					bank.ec.getShort('pindex',it) == index &&
					bank.ec.getFloat('lu',it).with{it < max_u && it > min_u} &&
					bank.ec.getFloat('lv',it).with{it < max_v && it > min_v} &&
					bank.ec.getFloat('lw',it).with{it < max_w && it > min_w} 
	}	
    }

    def passElectronEIEOCut= {bank, index ->
	//treat pcal as 'first layer' and ecal as second. cut on the edep in first layer
	return (0..<bank.ec.rows()).any{(bank.ec.getByte('detector',it) == DetectorType.ECAL.getDetectorId() &&
					 bank.ec.getInt('layer',it) == 1 &&
					 bank.ec.getFloat('energy',it) > min_pcal_dep &&
					 bank.ec.getShort('pindex',it) == index)
	}
    }

    //detector layer r1-12, r2-24, r3-36
    //rotate hit position based on sector
    def rotateDCHitPosition(hit,sec){
	def ang = Math.toRadians(sec*sect_angle_coverage)
	def x1_rot = hit.get(1) * Math.sin(ang) + hit.get(0) * Math.cos(ang)
	def y1_rot = hit.get(1) * Math.cos(ang) - hit.get(0) * Math.sin(ang)
	return [x1_rot,y1_rot]
    }
    
    //define left right 
    def borderDCHitPosition(y_rot,height){
	def slope = 1/Math.tan(Math.toRadians(0.5*sect_angle_coverage))
	def left  = (height - slope * y_rot)
	def right = (height + slope * y_rot)
	return [left, right]
    }

    def passElectronDCR1= {bank, index ->	
	def sec = (0..<bank.trck.rows()).findResult{(bank.trck.getShort('pindex',it) == index && 
						     bank.trck.getByte('detector',it) == DetectorType.DC.getDetectorId()) ? bank.trck.getByte('sector',it) : null }

	def hit_pos = (0..<bank.traj.rows()).findResult{(bank.traj.getShort('pindex',it) == index &&
							 bank.traj.getByte('detector',it) == DetectorType.DC.getDetectorId() &&
							 bank.traj.getByte('layer',it) == 12 )
							? 'xy'.collect{ axis -> bank.traj.getFloat(axis,it) } : null }

	if( sec == null || hit_pos == null ) return false	
		
	//get the sector for the track as defined in the REC::Track bank
	//x - 0 	//y - 1
	def hit_rotate = rotateDCHitPosition(hit_pos,sec-1)
	def left_right = borderDCHitPosition(hit_rotate.get(1),heightR1)	

	return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1) )  // && x1_rot**2 > radius2_DCr1){ return true }
	
    }

    def passElectronDCR2= {bank, index ->	
	def sec = (0..<bank.trck.rows()).findResult{(bank.trck.getShort('pindex',it) == index && 
					      bank.trck.getByte('detector',it) == DetectorType.DC.getDetectorId()) ? bank.trck.getByte('sector',it) : null }

	def hit_pos = (0..<bank.traj.rows()).findResult{(bank.traj.getShort('pindex',it) == index &&
							 bank.traj.getByte('detector',it) == DetectorType.DC.getDetectorId() &&
							 bank.traj.getByte('layer',it) == 24 )
							? 'xy'.collect{ axis -> bank.traj.getFloat(axis,it) } : null }

	if( sec == null || hit_pos == null ) return false	
		
	//get the sector for the track as defined in the REC::Track bank
	//x - 0 	//y - 1
	def hit_rotate = rotateDCHitPosition(hit_pos,sec-1)
	def left_right = borderDCHitPosition(hit_rotate.get(1),heightR2)	

	return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1) )  // && x1_rot**2 > radius2_DCr1){ return true }
	
    }

    def passElectronDCR3= {bank, index ->	
	def sec = (0..<bank.trck.rows()).findResult{(bank.trck.getShort('pindex',it) == index && 
					      bank.trck.getByte('detector',it) == DetectorType.DC.getDetectorId()) ? bank.trck.getByte('sector',it) : null }

	def hit_pos = (0..<bank.traj.rows()).findResult{(bank.traj.getShort('pindex',it) == index &&
							 bank.traj.getByte('detector',it) == DetectorType.DC.getDetectorId() &&
							 bank.traj.getByte('layer',it) == 36 )
							? 'xy'.collect{ axis -> bank.traj.getFloat(axis,it) } : null }

	if( sec == null || hit_pos == null ) return false	
	
	//get the sector for the track as defined in the REC::Track bank
	//x - 0 	//y - 1
	def hit_rotate = rotateDCHitPosition(hit_pos,sec-1)
	def left_right = borderDCHitPosition(hit_rotate.get(1),heightR3)		
	return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1) )  // && x1_rot**2 > radius2_DCr1){ return true }
	
    }


    // using the static def approach in the main code we do not need to instantiate the electron class -> use Electron.isElectronEBPIDCut
    //static def isElectronEBPIDCut(bank, iele){/
    //return partb.getInt('pid',iele) == 11 	
    //}

}
