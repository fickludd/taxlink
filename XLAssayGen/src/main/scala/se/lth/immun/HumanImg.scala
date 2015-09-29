package se.lth.immun

import java.io.File
import java.awt.Color
import java.awt.Graphics2D
import java.awt.FontMetrics
import java.awt.image.BufferedImage
import javax.imageio.ImageIO
import java.io.IOException

import se.lth.immun.xlink.XLink
import se.lth.immun.chem._
import se.lth.immun.unimod._

object HumanImg {

	import XLAssayGen._
	
	val w = 1000
	val h = 400
	
	def report(base:File, xlAssays:Seq[XLAssay]) = {
		for (xl <- xlAssays)
			reportXLAssay(base, xl)
	}
	
	def reportXLAssay(base:File, xl:XLAssay) = {
		
		val image = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB)
		val g = image.createGraphics
		
		val isCross = xl.xlink.positions.exists(l => l.pepPos1.pep != l.pepPos2.pep)
		if (isCross) {
			drawCross(g, xl)
		} else {
			drawLoop(g, xl)
		}
		
		g.dispose
		try { 
		    ImageIO.write(image, "png", new File(base.getAbsolutePath + xl.xlinkStr + ".png")) 
		} catch {
			case ioe:IOException =>
		    	ioe.printStackTrace
		}
	}
	
	
	def drawCross(g:Graphics2D, xl:XLAssay) = {
		
		import XLink._
		
		val fm = g.getFontMetrics
		val top = xl.xlink.pep1
		val bot = xl.xlink.pep2.get
		val topMap = new PeptideMap(top, true, fm)
		val botMap = new PeptideMap(bot, false, fm)
		
		val positions = xl.xlink.positions.flatMap(link => 
				List(link.pepPos1, link.pepPos2)
			).collect { 
				case kp:KnownPosition => kp 
			}
		
		val topPos = positions.filter(_.pep == top).map(_.pos)
		val botPos = positions.filter(_.pep == bot).map(_.pos)
				
		g.drawString(topMap.coreStr, topMap.x, topMap.y)
		g.drawString(botMap.coreStr, botMap.x, botMap.y)
		
		for (l <- xl.xlink.positions) {
			(l.pepPos1, l.pepPos2) match {
				case (KnownPosition(pos1, pep1), KnownPosition(pos2, pep2)) =>
					val m1 = if (pep1 == top) topMap else botMap
					val m2 = if (pep2 == top) topMap else botMap
					g.drawLine(m1(pos1) - 3, m1.yanchor, m2(pos2) - 3, m2.yanchor)
				case _ => {}
			}
		}
		
		for (t <- xl.transitions) {
			val m = if (t.peptide.top) topMap else botMap
			t.fragmentType match {
				case 'y' =>
					g.setColor(Color.RED)
					g.drawLine(m(t.pos), m.y - 20, m(t.pos), m.y)
				case 'b' =>
					g.setColor(Color.BLUE)
					g.drawLine(m(t.pos), m.y - 10, m(t.pos), m.y + 10)
				case _ => println("Unsupported img generation for fragment type "+t.fragmentType)
			}
		}
	}
	
	case class AA(stdAA:StandardAminoAcid, mod:String)
	def pepToAAs(pep:Peptide) = {
		def makeModStr(mod:IMolecule) =
			if (mod.toString.startsWith("UniMod:")) mod.toString
			else "%.2f".format(mod.monoisotopicMass)
		
		for (aa <- pep.aminoAcids) yield
			aa match {
				case saa:StandardAminoAcid =>
					AA(saa, "")
					
				case maa:ModifiedAminoAcid =>
					maa.aa match {
						case saa:StandardAminoAcid =>
							AA(saa, makeModStr(maa.modification))
						case _ =>
							throw new Exception("No support for nested modifications: "+maa)
					}
			}
	}
	
	
	class PeptideMap(val pep:Peptide, up:Boolean, fm:FontMetrics) {
		val y = if (up) h / 3 else (h*2) / 3
		val aas = pepToAAs(pep)
		val coreStr = aas.map(_.stdAA.letter).mkString
		val coreWidth = fm.charsWidth(coreStr.toArray, 0, coreStr.length)
		val x = w / 2 - coreWidth / 2
		
		def apply(pos:Int) =
			x + fm.charsWidth(coreStr.toArray, 0, pos)
			
		def yanchor =
			if (up) y + 10 else y - 20
	}
	
	
	def drawLoop(g:Graphics2D, xl:XLAssay) = {
		
	}
}