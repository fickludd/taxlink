package se.lth.immun

import java.io.File
import java.io.FileWriter
import java.io.BufferedWriter

import se.lth.immun.xml.XmlWriter
import se.lth.immun.chem._
import se.lth.immun.traml.ghost._

object XLTraml {

	import XLAssayGen._
	
	def write(f:File, xlAssays:Seq[XLAssay]) = {
		val gt = new GhostTraML
		
		for (xl <- xlAssays) {
			if (!gt.compounds.contains(xl.xlinkStr)) {
				val gc = new GhostCompound
				gc.id = xl.xlinkStr
				gc.mass = xl.mz * xl.z - Constants.PROTON_WEIGHT
				gt.compounds(gc.id) = gc
			}
			gt.compounds(xl.xlinkStr).preferredCharges = xl.z +: gt.compounds(xl.xlinkStr).preferredCharges
			
			for (t <- xl.targets) {
				val g = new GhostTarget
				g.compoundRef = xl.xlinkStr
				g.id = xl.xlinkStr + " %s isotope=%s naturally=%.1f%%".format(t.label, t.order, t.occurence*100) 
				g.intensity = t.occurence
				g.q1 = t.mz
				g.q1z = xl.z
				
				gt.includes += g
			}
			
			for (t <- xl.transitions) {
				val g = new GhostTransition
				g.compoundRef = xl.xlinkStr
				g.id = xl.xlinkStr + " %s%d in %s".format(t.fragmentType, t.ordinal, t.peptide)
				g.q1 = t.precMz
				g.q1z = xl.z
				g.q3 = t.fragMz
				g.q3z = 1
				g.ions += "%s%d".format(t.fragmentType, t.ordinal)
				
				gt.transitions += g
			}
		}
		
		gt.write(new XmlWriter(new BufferedWriter(new FileWriter(f))))
		
	}
}