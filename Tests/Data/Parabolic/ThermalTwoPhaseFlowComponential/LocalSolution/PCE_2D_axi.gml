<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.opengeosys.org/images/xsd/OpenGeoSysGLI.xsd" xmlns:ogs="http://www.opengeosys.org">
    <name>geometry</name>
    <points>
        <point id="0" x="0" y="0" z="0"/>
        <point id="1" x="1.5" y="0" z="0"/>
        <point id="2" x="1.5" y="2" z="0"/>
        <point id="3" x="0" y="2" z="0"/>
        <point id="4" x="0" y="0.5" z="0"/>
    </points>
    <polylines>
        <polyline id="0" name="bottom">
			<pnt>1</pnt>
			<pnt>0</pnt>
        </polyline>
		<polyline id="1" name="lateral">
			<pnt>2</pnt>
			<pnt>1</pnt>
        </polyline>	
		<polyline id="2" name="top">
			<pnt>2</pnt>
			<pnt>3</pnt>
        </polyline>
		<polyline id="3" name="heat_source">
			<pnt>4</pnt>
			<pnt>3</pnt>
        </polyline>		
    </polylines>
</OpenGeoSysGLI>
