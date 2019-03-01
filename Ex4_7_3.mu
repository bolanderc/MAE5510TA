{
	"MachUp": {
		"version": 4
	},
	"metadata": {
		"version": 4.4,
		"type": "Object",
		"generator": "Object3D.toJSON"
	},
	"geometries": [
		{
			"uuid": "4F63B009-1A45-425F-9B81-C4FA7DF1F00F",
			"type": "SphereBufferGeometry",
			"radius": 0.1,
			"widthSegments": 32,
			"heightSegments": 16,
			"phiStart": 0,
			"phiLength": 6.283185307179586,
			"thetaStart": 0,
			"thetaLength": 3.141592653589793
		},
		{
			"uuid": "77B2BA79-0A2A-48C1-B156-2B65BB937A10",
			"type": "SphereBufferGeometry",
			"radius": 0.1,
			"widthSegments": 32,
			"heightSegments": 16,
			"phiStart": 0,
			"phiLength": 6.283185307179586,
			"thetaStart": 0,
			"thetaLength": 3.141592653589793
		},
		{
			"uuid": "7E97483C-C2FE-4273-AA1C-E7FB918B05DE",
			"type": "WingGeometry",
			"is_main": true,
			"side": "both",
			"span": 16.5,
			"sweep": 0,
			"dihedral": 0,
			"mount": 1.8,
			"washout": 0,
			"root_chord": 7.7922,
			"tip_chord": 3.1169,
			"root_airfoil": {
				"NACA 2412": {
					"properties": {
						"type": "linear",
						"alpha_L0": -0.03838,
						"CL_alpha": 4.44,
						"Cm_L0": -0.0527,
						"Cm_alpha": -0.08,
						"CD0": 0.008,
						"CD0_L": 0,
						"CD0_L2": 0,
						"CL_max": 1.4
					}
				}
			},
			"tip_airfoil": {
				"NACA 2412": {
					"properties": {
						"type": "linear",
						"alpha_L0": -0.03838,
						"CL_alpha": 4.44,
						"Cm_L0": -0.0527,
						"Cm_alpha": -0.08,
						"CD0": 0.008,
						"CD0_L": 0,
						"CD0_L2": 0,
						"CL_max": 1.4
					}
				}
			},
			"nSeg": 40,
			"nAFseg": 50,
			"left_start": {
				"x": 0,
				"y": 0,
				"z": 0
			},
			"right_start": {
				"x": 0,
				"y": 0,
				"z": 0
			},
			"dy": 0,
			"control": {
				"has_control_surface": false,
				"span_root": 0.2,
				"span_tip": 0.8,
				"chord_root": 0.2,
				"chord_tip": 0.2,
				"is_sealed": 1,
				"mix": {
					"elevator": 1,
					"rudder": 0,
					"aileron": 0,
					"flap": 0
				}
			},
			"same_as_root": true
		},
		{
			"uuid": "7614D1BB-8238-4CFD-9B41-2BC513513ACD",
			"type": "WingGeometry",
			"is_main": false,
			"side": "both",
			"span": 7.379,
			"sweep": 0,
			"dihedral": 0,
			"mount": 4,
			"washout": 0,
			"root_chord": 3.4848,
			"tip_chord": 1.3939,
			"root_airfoil": {
				"NACA 2412": {
					"properties": {
						"type": "linear",
						"alpha_L0": -0.03838,
						"CL_alpha": 4.44,
						"Cm_L0": -0.0527,
						"Cm_alpha": -0.08,
						"CD0": 0.008,
						"CD0_L": 0,
						"CD0_L2": 0,
						"CL_max": 1.4
					}
				}
			},
			"tip_airfoil": {
				"NACA 2412": {
					"properties": {
						"type": "linear",
						"alpha_L0": -0.03838,
						"CL_alpha": 4.44,
						"Cm_L0": -0.0527,
						"Cm_alpha": -0.08,
						"CD0": 0.008,
						"CD0_L": 0,
						"CD0_L2": 0,
						"CL_max": 1.4
					}
				}
			},
			"nSeg": 40,
			"nAFseg": 50,
			"left_start": {
				"x": 0,
				"y": 0,
				"z": 0
			},
			"right_start": {
				"x": 0,
				"y": 0,
				"z": 0
			},
			"dy": 0,
			"control": {
				"has_control_surface": false,
				"span_root": 0.2,
				"span_tip": 0.8,
				"chord_root": 0.2,
				"chord_tip": 0.2,
				"is_sealed": 1,
				"mix": {
					"elevator": 1,
					"rudder": 0,
					"aileron": 0,
					"flap": 0
				}
			},
			"same_as_root": true
		}],
	"materials": [
		{
			"uuid": "E5893532-76A6-434E-948A-FD9B75C8F0D3",
			"type": "MeshStandardMaterial",
			"color": 16711680,
			"roughness": 0.5,
			"metalness": 0.5,
			"emissive": 16711680,
			"side": 2,
			"depthFunc": 3,
			"depthTest": true,
			"depthWrite": true,
			"skinning": false,
			"morphTargets": false
		},
		{
			"uuid": "CE9B2CE8-0573-4622-9A5B-A6D06AF21AEC",
			"type": "MeshStandardMaterial",
			"color": 6684927,
			"roughness": 0.5,
			"metalness": 0.5,
			"emissive": 6684927,
			"side": 2,
			"depthFunc": 3,
			"depthTest": true,
			"depthWrite": true,
			"skinning": false,
			"morphTargets": false
		},
		{
			"uuid": "D8EDC724-8A15-4CED-8051-0ACFC5C2DB30",
			"type": "MeshPhongMaterial",
			"color": 16777215,
			"emissive": 0,
			"specular": 1118481,
			"shininess": 30,
			"side": 2,
			"depthFunc": 3,
			"depthTest": true,
			"depthWrite": true,
			"skinning": false,
			"morphTargets": false
		},
		{
			"uuid": "5D43C091-1E1A-431B-8491-679A7684A9F0",
			"type": "MeshPhongMaterial",
			"color": 16777215,
			"emissive": 0,
			"specular": 1118481,
			"shininess": 30,
			"side": 2,
			"depthFunc": 3,
			"depthTest": true,
			"depthWrite": true,
			"skinning": false,
			"morphTargets": false
		}],
	"object": {
		"uuid": "84EB89B7-7044-4F3B-9BEE-E3BE072A11B0",
		"type": "Origin",
		"name": "MyAirplane",
		"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
		"children": [
			{
				"uuid": "98326C22-F088-48FF-A79A-639508D2CBD3",
				"type": "Mesh",
				"name": "Center of Gravity",
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
				"geometry": "4F63B009-1A45-425F-9B81-C4FA7DF1F00F",
				"material": "E5893532-76A6-434E-948A-FD9B75C8F0D3"
			},
			{
				"uuid": "E113E0CF-1E3B-4A7E-8FF6-85A24D64F658",
				"type": "Mesh",
				"name": "Aerodynamic Center",
				"visible": false,
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
				"geometry": "77B2BA79-0A2A-48C1-B156-2B65BB937A10",
				"material": "CE9B2CE8-0573-4622-9A5B-A6D06AF21AEC"
			},
			{
				"uuid": "F4900440-A077-4B1B-BE3F-CDC1D7D4D6C6",
				"type": "PointLight",
				"name": "PointLight",
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,10,10,-20,1],
				"color": 16777215,
				"intensity": 1,
				"distance": 0,
				"decay": 1
			},
			{
				"uuid": "1D3BD661-3FC8-4F30-B26D-7CE5E5173D50",
				"type": "Mesh",
				"name": "Main",
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,-3,0,-4,1],
				"geometry": "7E97483C-C2FE-4273-AA1C-E7FB918B05DE",
				"material": "D8EDC724-8A15-4CED-8051-0ACFC5C2DB30"
			},
			{
				"uuid": "8238A1E7-4451-48AA-B60C-91FA393BA2E0",
				"type": "Mesh",
				"name": "Canard",
				"matrix": [1,0,0,0,0,1,0,0,0,0,1,0,12,0,0,1],
				"geometry": "7614D1BB-8238-4CFD-9B41-2BC513513ACD",
				"material": "5D43C091-1E1A-431B-8491-679A7684A9F0"
			}],
		"background": 11184810
	}
}