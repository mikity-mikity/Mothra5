using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mikity.ghComponents
{
    public partial class Mothra4 : Grasshopper.Kernel.GH_Component
    {
        public override void BakeGeometry(Rhino.RhinoDoc doc, Rhino.DocObjects.ObjectAttributes att, List<Guid> obj_ids)
        {
            Rhino.Geometry.Transform zDown_airy = Rhino.Geometry.Transform.Translation(0, 0, 2d);
            Rhino.DocObjects.ObjectAttributes a2 = att.Duplicate();
            //a2.LayerIndex = 2;
            Rhino.DocObjects.ObjectAttributes a3 = att.Duplicate();
            //a3.LayerIndex = 3;
            Rhino.DocObjects.ObjectAttributes a4 = att.Duplicate();
            //a4.LayerIndex = 4;
            Rhino.DocObjects.ObjectAttributes a5 = att.Duplicate();
            //a5.LayerIndex = 5;
            Rhino.DocObjects.ObjectAttributes a6 = att.Duplicate();
            //a6.LayerIndex = 6;
            Rhino.DocObjects.ObjectAttributes a7 = att.Duplicate();
            //a7.LayerIndex = 7;
            foreach (var leaf in listLeaf)
            {
                var airySrf=leaf.airySrf.Duplicate() as Rhino.Geometry.NurbsSurface;
                airySrf.Transform(zScale);
                airySrf.Transform(zDown_airy);
                Guid id = doc.Objects.AddSurface(airySrf, a2);
                obj_ids.Add(id);
                
                var srf = leaf.shellSrf.Duplicate() as Rhino.Geometry.NurbsSurface;
                //srf.Transform(zDown_eq);
                id = doc.Objects.AddSurface(srf, a3);
                obj_ids.Add(id);
            }
            foreach (var branch in listBranch)
            {
                //var airyCrv=branch.airyCrv.Duplicate() as Rhino.Geometry.NurbsCurve;
                //airyCrv.Transform(zDown_airy);
                //Guid id = doc.Objects.AddCurve(airyCrv, a2);
                //obj_ids.Add(id);
                if (branch.branchType == branch.type.kink || branch.branchType == branch.type.reinforce||branch.branchType == branch.type.open)
                {
                    var crv = branch.shellCrv.Duplicate() as Rhino.Geometry.NurbsCurve;
                  //  crv.Transform(zDown_eq);
                    Guid id = doc.Objects.AddCurve(crv, a7);
                    obj_ids.Add(id);
                }
            }
            /*if (crossMagenta != null)
            {
                foreach (var line in crossMagenta)
                {
                    Guid id = doc.Objects.AddLine(line, a4);
                    obj_ids.Add(id);
                }
            }
            if (listBranch != null)
            {
                foreach (var branch in listBranch)
                {
                    if (branch.branchType == branch.type.fix)
                    {
                        if (branch.tuples != null)
                        {
                            foreach (var tup in branch.tuples)
                            {
                                var circle = new Rhino.Geometry.Circle(new Rhino.Geometry.Point3d(tup.x, tup.y, tup.z), 0.5);
                                circle.Transform(zDown);
                                Guid id = doc.Objects.AddCircle(circle, a6);
                                obj_ids.Add(id);
                                circle = new Rhino.Geometry.Circle(new Rhino.Geometry.Point3d(tup.x, tup.y, branch.slice2.height), 0.5);
                                circle.Transform(zDown_eq);
                                id = doc.Objects.AddCircle(circle, a6);
                                obj_ids.Add(id);
                            }
                        }
                    }
                    if (branch.branchType == branch.type.kink || branch.branchType == branch.type.reinforce||branch.branchType==branch.type.open)
                    {
                        if (branch.tuples != null)
                        {
                            foreach (var tup in branch.tuples)
                            {
                                var D = tup.SPK[0, 0]/2d;
                                if (D > 0)
                                {
                                    var line = new Rhino.Geometry.Line(new Rhino.Geometry.Point3d(tup.x, tup.y, tup.z), new Rhino.Geometry.Point3d(tup.x, tup.y, tup.z + D));
                                    line.Transform(zDown);
                                    Guid id = doc.Objects.AddLine(line, a5);
                                    obj_ids.Add(id);
                                }
                            }
                        }
                    }
                    
                }
            }*/
        }
    }
}
