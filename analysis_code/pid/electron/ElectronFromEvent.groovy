package pid.electron

import org.jlab.detector.base.DetectorType
import event.Event

class ElectronFromEvent {

    def ebPID = 11
    def min_nphe = 2

    //vertex wide enough for all sectors
    def min_vz = -12
    def max_vz = 9

    def min_ecal_dep = 0.06

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

    def passElectronStatus = { event, index ->
        return (event.status[index] < 0)
    }

    def passElectronEBPIDCut = { event, index ->
		return (event.pid[index] == ebPID)
    }

    def passElectronChargeCut = { event, index ->
		return (event.charge[index] < 0)
    }

    def passElectronNpheCut = { event, index ->
		return (event.cherenkov_status.contains(index)) ? event.nphe[index] > min_nphe : false
    }

    def passElectronVertexCut = { event, index ->
        return (event.vz[index].with{ it < max_vz && it > min_vz })
    }

    def passElectronPCALFiducialCut = { event, index ->
		if (event.pcal_status.contains(index)){
			return (event.pcal_u[index].with{it < max_u && it > min_u } &&
					event.pcal_v[index].with{it < max_v && it > min_v } &&
					event.pcal_w[index].with{it < max_w && it > min_w })
		}
		return false
    }

    def passElectronEIEOCut = { event, index ->
		return (event.ecal_inner_status.contains(index)) ? event.ecal_inner_energy[index] > min_ecal_dep : false
	}

    //detector layer r1-12, r2-24, r3-36
    //rotate hit position based on sector
    def rotateDCHitPosition(hit, sec) {
        def ang = Math.toRadians(sec * sect_angle_coverage)
        def x1_rot = hit.get(1) * Math.sin(ang) + hit.get(0) * Math.cos(ang)
        def y1_rot = hit.get(1) * Math.cos(ang) - hit.get(0) * Math.sin(ang)
        return [x1_rot, y1_rot]
    }

    //define left right 
    def borderDCHitPosition(y_rot, height) {
        def slope = 1 / Math.tan(Math.toRadians(0.5 * sect_angle_coverage))
        def left = (height - slope * y_rot)
        def right = (height + slope * y_rot)
        return [left, right]
    }

	def passElectronDCR1 = { event, index ->
		if (event.dc1_status.contains(index)){
			def sec = event.dc_sector[index]
			def hit = event.dc1.get(index).find{ hit -> hit.layer == 12}
            if (hit){
                def hit_rotate = rotateDCHitPosition([hit.x, hit.y], sec-1)
                def left_right = borderDCHitPosition(hit_rotate.get(1), heightR1)
                return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1))
            } else {
                return false
            }
		}
		return false
	}

	def passElectronDCR2 = { event, index ->
		if (event.dc2_status.contains(index)){
			def sec = event.dc_sector[index]
			def hit = event.dc2.get(index).find{ hit -> hit.layer == 24}
            if (hit){
                def hit_rotate = rotateDCHitPosition([hit.x, hit.y], sec-1)
                def left_right = borderDCHitPosition(hit_rotate.get(1), heightR2)
                return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1))
            } else {
                return false
            }
		}
		return false
	}

	def passElectronDCR3 = { event, index ->
		if (event.dc3_status.contains(index)){
			def sec = event.dc_sector[index]
			def hit = event.dc3.get(index).find{ hit -> hit.layer == 36}
			if (hit) {
                def hit_rotate = rotateDCHitPosition([hit.x, hit.y], sec-1)
                def left_right = borderDCHitPosition(hit_rotate.get(1), heightR3)

                // from lower line: && x1_rot**2 > radius2_DCr1){ return true }
                return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1))
            } else {
                return false
            }
		}
		return false
	}
}
