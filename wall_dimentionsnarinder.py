
>>> import ifcopenshell
>>> ff = ifcopenshell.open("C:\Users\narindersingh\Desktop\Narinder\sdfghj.ifc")
>>> wall = ff.by_type("ifcwall")
>>> wall = ff.by_type("ifcwall")
>>> print wall
>>> wlen = wall.Representation.Representations[0].Items[0].SweptArea.XDim
>>> print wlen
>>> whgt = wall.Representation.Representations[0].Items[0].SweptArea.YDim
>>> wht = wall.Representation.Representations[0].Items[0].Depth
>>> print wht

 
