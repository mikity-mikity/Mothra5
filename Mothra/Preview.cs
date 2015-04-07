using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mikity.ghComponents
{
    public partial class Mothra4 : Grasshopper.Kernel.GH_Component
    {
        public Rhino.Geometry.Transform zDown_eq = Rhino.Geometry.Transform.Translation(0, 0, 25d);
        public Rhino.Geometry.Transform zDown = Rhino.Geometry.Transform.Translation(0, 0, 15d);
        public Rhino.Geometry.Transform zScale = Rhino.Geometry.Transform.Scale(Rhino.Geometry.Plane.WorldXY, 1, 1, 1d);
        public override void DrawViewportWires(Grasshopper.Kernel.IGH_PreviewArgs args)
        {
            if (Hidden)
            {
                return;
            }
            if (listArrow != null)
            {
                args.Display.DrawLines(listArrow, System.Drawing.Color.Red);
            }
            //targetSrf
            if (listPnt != null)
            {
                args.Display.DrawPoints(targetSrf, Rhino.Display.PointStyle.Simple, 1, System.Drawing.Color.White);
            }
            //eigenvectors
            if (crossCyan != null)
            {
                args.Display.DrawLines(crossCyan, System.Drawing.Color.Cyan);
            }
            if (crossMagenta != null)
            {
                args.Display.DrawLines(crossMagenta, System.Drawing.Color.Magenta);
            }
            if (listError != null)
            {
                foreach (var error in listError)
                {
                    args.Display.DrawCurve(error, System.Drawing.Color.Red, 10);
                }
            }
            /*if (listSlice != null)
            {
                foreach (var slice in listSlice.Values)
                {
                    var pl = slice.pl;
                    args.Display.DrawPolygon(new Rhino.Geometry.Point3d[] { pl.PointAt(-2, -2), pl.PointAt(-2, 2), pl.PointAt(2, 2), pl.PointAt(2, -2) }, System.Drawing.Color.Azure, true);
                }
            }*/
            foreach (var branch in listBranch)
            {
                switch (branch.branchType)
                {
                    case Mothra4.branch.type.fix:
                        args.Display.DrawCurve(branch.crv, System.Drawing.Color.Orange, 3);
                        break;
                    case Mothra4.branch.type.reinforce:
                        args.Display.DrawCurve(branch.crv, System.Drawing.Color.Cyan, 3);
                        break;
                    case Mothra4.branch.type.kink:
                        args.Display.DrawCurve(branch.crv, System.Drawing.Color.Purple, 3);
                        break;
                    case Mothra4.branch.type.open:
                        args.Display.DrawCurve(branch.crv, System.Drawing.Color.Green, 3);
                        break;
                }
            }
            foreach (var leaf in listLeaf)
            {
                if (leaf.airySrf != null)
                {
                    var srf = leaf.airySrf.Duplicate() as Rhino.Geometry.NurbsSurface;
                    srf.Transform(zScale);
                    args.Display.DrawSurface(srf, System.Drawing.Color.Brown, 3);
                }
            }
            foreach (var leaf in listLeaf)
            {
                if (leaf.shellSrf != null)
                {
                    var srf = leaf.shellSrf.Duplicate() as Rhino.Geometry.NurbsSurface;
                    srf.Transform(zDown_eq);
                    args.Display.DrawSurface(srf, System.Drawing.Color.Brown, 3);
                }
            }
            /*foreach (var branch in listBranch)
            {
                if (branch.airyCrv != null)
                {
                    var crv = branch.airyCrv.Duplicate() as Rhino.Geometry.NurbsCurve;
                    crv.Transform(zScale);
                    args.Display.DrawCurve(crv, System.Drawing.Color.SeaGreen, 4);
                }
            }*/
            foreach (var branch in listBranch)
            {
                if (branch.branchType == branch.type.fix)
                {
                    double x = 0, y = 0, z = branch.slice2.height;
                    for (int i = 0; i < branch.N; i++)
                    {
                        x += branch.crv.Points[i].Location.X;
                        y += branch.crv.Points[i].Location.Y;
                        z += branch.crv.Points[i].Location.Z;
                    }
                    x /= branch.N;
                    y /= branch.N;
                    args.Display.DrawLine(new Rhino.Geometry.Point3d(x - 1d, y - 1d, z), new Rhino.Geometry.Point3d(x + 1d, y - 1d, z), System.Drawing.Color.Blue);
                    args.Display.DrawLine(new Rhino.Geometry.Point3d(x - 1d, y - 1d, z), new Rhino.Geometry.Point3d(x - 1d, y + 1d, z), System.Drawing.Color.Blue);
                    args.Display.DrawLine(new Rhino.Geometry.Point3d(x + 1d, y + 1d, z), new Rhino.Geometry.Point3d(x + 1d, y - 1d, z), System.Drawing.Color.Blue);
                    args.Display.DrawLine(new Rhino.Geometry.Point3d(x + 1d, y + 1d, z), new Rhino.Geometry.Point3d(x - 1d, y + 1d, z), System.Drawing.Color.Blue);
                }
            }
            foreach (var branch in listBranch)
            {
                if (branch.shellCrv != null)
                {
                    if (branch.branchType == branch.type.kink || branch.branchType == branch.type.reinforce || branch.branchType == branch.type.open)
                    {

                        var crv = branch.shellCrv.Duplicate() as Rhino.Geometry.NurbsCurve;
                        crv.Transform(zDown_eq);
                        args.Display.DrawCurve(crv, System.Drawing.Color.SeaGreen, 3);
                    }
                }
            }
            //find max value of reinforcement
            if (listBranch != null)
            {
                foreach (var branch in listBranch)
                {
                    if (branch.tuples != null)
                    {
                        if (branch.branchType == branch.type.fix)
                        {
                            foreach (var tup in branch.tuples)
                            {
                                var circle = new Rhino.Geometry.Circle(new Rhino.Geometry.Point3d(tup.x, tup.y, tup.z), 0.5);
                                circle.Transform(zDown);
                                args.Display.DrawCircle(circle, System.Drawing.Color.Yellow, 2);
                                circle = new Rhino.Geometry.Circle(new Rhino.Geometry.Point3d(tup.x, tup.y, branch.slice2.height), 0.5);
                                circle.Transform(zDown_eq);
                                args.Display.DrawCircle(circle, System.Drawing.Color.Yellow, 2);
                            }
                        }
                        else
                        {
                            foreach (var tup in branch.tuples)
                            {
                                var D = tup.SPK[0, 0]/2d;
                                if (D > 0)
                                {
                                    var line = new Rhino.Geometry.Line(new Rhino.Geometry.Point3d(tup.x, tup.y, tup.z), new Rhino.Geometry.Point3d(tup.x, tup.y, tup.z + D));
                                    line.Transform(zDown);
                                    args.Display.DrawLine(line, System.Drawing.Color.Red, 2);
                                }
                                else
                                {
                                    var line = new Rhino.Geometry.Line(new Rhino.Geometry.Point3d(tup.x, tup.y, tup.z), new Rhino.Geometry.Point3d(tup.x, tup.y, tup.z + D));
                                    line.Transform(zDown);
                                    args.Display.DrawLine(line, System.Drawing.Color.Blue, 2);
                                }
                            }
                        }
                    }
                }
            }
    
            /*
            if (a != null)
            {
                args.Display.DrawPoints(a, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Blue);

            }*/
            if (a2 != null)
            {
                args.Display.DrawPoints(a2, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Blue);

            }
        }
    }
}
