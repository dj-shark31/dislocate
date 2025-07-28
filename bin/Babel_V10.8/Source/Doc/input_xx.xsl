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
						<a href="./input_prepareDrag.html">prepareDrag</a>,
						<a href="./structures.html">structures</a>
					</font>
				</div>
				<!--Title of the web page-->
				<h1 id="headTitle">
					<xsl:value-of select="manual/title"/>
				</h1>
				<!--List of the namelists-->
				<xsl:if test="count(manual/namelist) &gt; 1">
				<i>Namelists: </i>
				<xsl:for-each select="manual/namelist">
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
				<!--Namelists-->
				<xsl:apply-templates select="manual/namelist"/>
				<!--Footnote-->
				<div style="text-align:right">
					<font size="2">
						Back to 
						<a href="./index.html">main page</a>,
						<a href="./input_babel.html">babel</a>,
						<a href="./input_displacement.html">displacement</a>,
						<a href="./input_patternDetect.html">patternDetect</a>,
						<a href="./input_prepareDrag.html">prepareDrag</a>,
						<a href="./structures.html">structures</a>
					</font>
				</div>
			</BODY>
		</HTML>
	</xsl:template>

	<!--Namelist template-->
	<xsl:template match="namelist">
		<!--Namelist name-->
		<h2>
			<xsl:attribute name="id"><xsl:value-of select="normalize-space(name)"/></xsl:attribute>
			Namelist <xsl:value-of select="name"/>
			<!--Link to the head of the page -->
			<a href="#headTitle">&#x2191;</a>
		</h2>
		<!--List of variables: -->
		<i>Variables: </i>
		<xsl:for-each select="var">
				<a>
				<xsl:attribute name="href">#<xsl:value-of select="normalize-space(name)"/></xsl:attribute>
				<xsl:value-of select="normalize-space(name)"/>
				</a>,
		</xsl:for-each>
		<!--Description of the namelist-->
		<div style="white-space: pre-line;">
		<xsl:value-of select="info"/>
		</div>
		<!--Description of all variables-->
		<ul>
			<xsl:apply-templates select="var"/>
		</ul>
	</xsl:template>

	<!--Var template-->
	<xsl:template match="var">
		<!--Variable name-->
		<b>
			<xsl:attribute name="id"><xsl:value-of select="normalize-space(name)"/></xsl:attribute>
			<xsl:value-of select="name"/> 
		</b>
		<!--Link to the parent namelist-->
		<a>
			<xsl:attribute name="href">#<xsl:value-of select="normalize-space(../name)"/></xsl:attribute>
			&#x2191;
		</a>
		<!--Variable description-->
		<br/>
		<ul>
			<!--Properties-->
			(<i>type</i>: <xsl:value-of select="type"/>, 
			<i>default</i>: <xsl:value-of select="default"/>,
			<i>mandatory</i>: <xsl:value-of select="mandatory"/>)
			<div style="white-space: pre-line;">
			<xsl:value-of select="info"/>
			</div>
			<!--Related variables-->
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
		</ul>
		<br/> <br/>
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

