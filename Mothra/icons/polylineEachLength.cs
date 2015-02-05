using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using mikity;
using mikity.LinearAlgebra;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using mikity.NumericalMethodHelper;
using mikity.NumericalMethodHelper.objects;
namespace mikity.ghComponents
{
    /// <summary>
    /// Construct a point array using isoparametric shape functions.
    /// </summary>
    public class polyline_each_length : Grasshopper.Kernel.GH_Component
    {
        /*protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //現在のコードを実行しているAssemblyを取得
                System.Reflection.Assembly myAssembly =
                    System.Reflection.Assembly.GetExecutingAssembly();

                System.IO.Stream st = myAssembly.GetManifestResourceStream("mikity.ghComponents.icons.icon34.bmp");
                //指定されたマニフェストリソースを読み込む
                System.Drawing.Bitmap bmp = new System.Drawing.Bitmap(st);
                return bmp;
            }
        }*/

        public polyline_each_length()
            : base("polyline->eachLength [Constraint Condition]", "polyline->eachLength", "polyline->eachLength", "Ricecooker", "Super market")
        {
        }
        protected override void RegisterInputParams(Grasshopper.Kernel.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Curve", "c", "Curve", Grasshopper.Kernel.GH_ParamAccess.item);
            pManager.AddNumberParameter("Length", "L", "Prescribed value for the length", Grasshopper.Kernel.GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(Grasshopper.Kernel.GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Particle System", "pS", "Particle System", Grasshopper.Kernel.GH_ParamAccess.item);
            pManager.Register_PointParam("Points", "P", "Particle System", Grasshopper.Kernel.GH_ParamAccess.list);
        }
        public override void DrawViewportWires(Grasshopper.Kernel.IGH_PreviewArgs args)
        {
            if (Hidden)
            {
                return;
            }
            if (this.Attributes.Selected)
            {
                args.Display.DrawLines(lGeometry, System.Drawing.Color.Magenta, 3);
            }
            else
            {
                args.Display.DrawLines(lGeometry, System.Drawing.Color.DarkMagenta, 3);
            }
            base.DrawViewportWires(args);
        }
        private const int _nNodes = 2;
        private const int _dim = 1;

        GH_particleSystem pS;
        int nNewNodes = 0;
        int nElements = 0;
        List<mikity.NumericalMethodHelper.objects.constrainVolumeObject> lCV = new List<mikity.NumericalMethodHelper.objects.constrainVolumeObject>();
        List<Rhino.Geometry.Line> lGeometry = new List<Rhino.Geometry.Line>();
        List<Rhino.Geometry.Point3d> newNodes = new List<Rhino.Geometry.Point3d>();
        protected override void SolveInstance(Grasshopper.Kernel.IGH_DataAccess DA)
        {
            double v = -1.0;
            DA.GetData(1, ref v);
            if (!FriedChiken.isInitialized)
            {
                Rhino.Geometry.Curve c = null;
                if (!DA.GetData(0, ref c)) return;
                if (c.IsPolyline())
                {
                    Rhino.Geometry.Polyline pl = null;
                    if (c.TryGetPolyline(out pl))
                    {
                        nNewNodes = pl.Count();
                        nElements = nNewNodes - 1;
                        newNodes.Clear();
                        newNodes.AddRange(pl);
                        lGeometry.Clear();
                        for (int i = 0; i < nElements; i++)
                        {
                            lGeometry.Add(new Rhino.Geometry.Line(newNodes[i], newNodes[i + 1]));
                        }


                        mikity.NumericalMethodHelper.particle[] particles = new mikity.NumericalMethodHelper.particle[nNewNodes];
                        for (int i = 0; i < nNewNodes; i++)
                        {
                            particles[i] = new mikity.NumericalMethodHelper.particle(newNodes[i][0], newNodes[i][1], newNodes[i][2]);
                        }
                        List<mikity.NumericalMethodHelper.elements.isoparametricElement> e = new List<mikity.NumericalMethodHelper.elements.isoparametricElement>();
                        for (int i = 0; i < nElements; i++)
                        {
                            e.Add(new mikity.NumericalMethodHelper.elements.isoparametricElement(i, i + 1));
                        }
                        lCV.Clear();
                        pS = new GH_particleSystem(particles);
                        for (int i = 0; i < e.Count; i++)
                        {
                            lCV.Add(new constrainVolumeObject(v / nElements));
                            lCV[i].addElement(e[i]);
                            pS.Value.addObject(lCV[i]);
                        }

                        lGeometry.Clear();
                        for (int i = 0; i < nElements; i++)
                        {
                            lGeometry.Add(new Rhino.Geometry.Line(pS.Value.particles[i, 0], pS.Value.particles[i, 1], pS.Value.particles[i, 2], pS.Value.particles[i + 1, 0], pS.Value.particles[i + 1, 1], pS.Value.particles[i + 1, 2]));
                        }
                    }
                }
                else
                {
                    AddRuntimeMessage(Grasshopper.Kernel.GH_RuntimeMessageLevel.Error, "Only polyline is accepted");
                    return;
                }

            }
            else
            {
                if (lCV != null && v > 0)
                {
                    for (int i = 0; i < lCV.Count; i++)
                    {
                        lCV[i].refVolume = v / nElements;
                    }
                }
            }
            DA.SetData(0, pS);
            DA.SetDataList(1, newNodes);
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("56ae1ad7-757f-4fa3-b795-d2f627e4d405"); }
        }

    }
}
