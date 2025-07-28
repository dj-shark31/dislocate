<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
>
<xsl:output method='html' xml:space="preserve" encoding="UTF-8"/>
	<xsl:template match="/">
		<HTML>
			<HEAD>
				<TITLE>
					<xsl:value-of select="manual/title"/>
				</TITLE>
			</HEAD>
			<BODY BGCOLOR="#FFFFFF">
				<!--Header-->
				<div style="text-align:right">
					<font size="2">
						Back to 
						<a href="./index.html">main page</a>,
						<a href="./input_babel.html">babel</a>,
						<a href="./input_displacement.html">displacement</a>,
						<a href="./input_patternInit.html">patternInit</a>, 
						<a href="./input_patternDetect.html">patternDetect</a>, 
						<a href="./input_prepareDrag.html">prepareDrag</a>
					</font>
				</div>
				<!--Title of the web page-->
				<h1 id="headTitle">
					<xsl:value-of select="manual/title"/>
				</h1>
				<!--List of the structure types-->
				<xsl:if test="count(manual/type) &gt; 1">
				<i>Types: </i>
				<xsl:for-each select="manual/type">
						<a>
						<xsl:attribute name="href">#<xsl:value-of select="normalize-space(name)"/></xsl:attribute>
						<xsl:value-of select="normalize-space(name)"/>
						</a>,
				</xsl:for-each>
				</xsl:if>
				<!--Introduction-->
				<div style="white-space: pre-line;">
				<xsl:value-of select="manual/intro"/>
				</div>
				<!--Types-->
				<xsl:apply-templates select="manual/type"/>
				<!--Footnote-->
				<div style="text-align:right">
					<font size="2">
						Back to 
						<a href="./index.html">main page</a>,
						<a href="./input_babel.html">babel</a>,
						<a href="./input_displacement.html">displacement</a>,
						<a href="./input_patternDetect.html">patternDetect</a>
					</font>
				</div>
			</BODY>
		</HTML>
	</xsl:template>

	<!--Type template-->
	<xsl:template match="type">
		<!--Type name-->
		<h2>
			<xsl:attribute name="id"><xsl:value-of select="normalize-space(name)"/></xsl:attribute>
			<xsl:value-of select="name"/> type
			<!--Link to the head of the page -->
			<a href="#headTitle">&#x2191;</a>
		</h2>
		(<i>keywords</i>:
		<xsl:value-of select="keywords"/>)
		<!--Description of the type-->
		<div style="white-space: pre-line;">
		<xsl:value-of select="info"/>
		</div>
		<!--Related info-->
		<xsl:for-each select="see_also[1]">
			<br/>
			<i>See also: </i>
			<xsl:value-of select="text"/>
			<xsl:apply-templates select="ilink"/>
			<xsl:apply-templates select="elink"/>
		</xsl:for-each>
		<xsl:for-each select="see_also[position()>1]">
			,
			<xsl:value-of select="text"/>
			<xsl:apply-templates select="ilink"/>
			<xsl:apply-templates select="elink"/>
		</xsl:for-each>
	</xsl:template>

	<!--Internal link template-->
	<xsl:template match="ilink">
				<a>
				<xsl:attribute name="href">#<xsl:value-of select="normalize-space(.)"/></xsl:attribute>
				<xsl:value-of select="normalize-space(.)"/>
				</a>
	</xsl:template>
	<!--External link template-->
	<xsl:template match="elink">
				<a>
				<xsl:attribute name="href"><xsl:value-of select="normalize-space(url)"/></xsl:attribute>
				<xsl:value-of select="normalize-space(text)"/>
				</a>
	</xsl:template>

</xsl:stylesheet>

