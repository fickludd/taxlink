package se.lth.immun

class Target(
		val mz:Double,
		val nXLinks:Int,
		val charge:Int,
		val id:String
) {
	override def toString = "%.5f +%d %d-linked".format(mz, charge, nXLinks)
}