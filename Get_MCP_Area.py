import arcpy
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=1

##lyr = arcpy.KMLToLayer_conversion(kml, wkdir)

wkdir = r"C:\Users\mevans\GIS\LPC"
ft = r"%s\Andrews_EE_Wells2_prj.shp"%(wkdir)
nlcd = r"%s\NLCD11_TX"%(wkdir)
temp = r"%s\Temp"%(wkdir)
table = arcpy.sa.ZonalStatisticsAsTable(ft, "FID", nlcd, r"%s\table"%(temp), "DATA", "MAJORITY")
mcp = arcpy.MinimumBoundingGeometry_management(ft, r"%s\temp.shp"%(wkdir),
                                                "CONVEX_HULL")
arcpy.AddField_management(mcp, "PERIM", "FLOAT")
arcpy.AddField_management(mcp, "AREA", "FLOAT")
expression1 = "{0}".format("!SHAPE.area@ACRES!")
expression2 = "{0}".format("!SHAPE.length@METERS!")
arcpy.CalculateField_management(mcp, "AREA", expression1, "PYTHON")
arcpy.CalculateField_management(mcp, "PERIM", expression2, "PYTHON")
arcpy.JoinField_management(ft, "FID", mcp, "ORIG_FID", ["PERIM", "AREA"]) 
arcpy.JoinField_management(ft, "FID", table, "FID", ["MAJORITY"])
##cursor = arcpy.UpdateCursor(ft)
##for row in cursor:
## print(row.getValue("FID"))
##
## geom = mcp.GetGeometryRef()
## area = geom.GetArea() * 0.000247105
## row.setValue("MCP", area)
## cursor.updateRow(row)

#del row, cursor


