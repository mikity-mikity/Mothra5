using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mikity.ghComponents
{
    public partial class Mothra4 : Grasshopper.Kernel.GH_Component
    {
        void Nurbs2x(Rhino.Geometry.NurbsSurface srf, double[,] _x)
        {
            for (int i = 0; i < srf.Points.CountV; i++)
            {
                for (int j = 0; j < srf.Points.CountU; j++)
                {
                    _x[(i * srf.Points.CountU + j), 0] = srf.Points.GetControlPoint(j, i).Location.X;
                    _x[(i * srf.Points.CountU + j), 1] = srf.Points.GetControlPoint(j, i).Location.Y;
                    _x[(i * srf.Points.CountU + j), 2] = srf.Points.GetControlPoint(j, i).Location.Z;
                }
            }
        }
        void Nurbs2x(Rhino.Geometry.NurbsCurve crv, double[,] _x)
        {
            for (int i = 0; i < crv.Points.Count; i++)
            {
                _x[i, 0] = crv.Points[i].Location.X;
                _x[i, 1] = crv.Points[i].Location.Y;
                _x[i, 2] = crv.Points[i].Location.Z;
            }
        }
        void createNurbsElements(branch branch)
        {
            double[] uKnot;

            int N = branch.N;
            int uDim = branch.crv.Order;
            int uDdim = branch.crv.Order - 1;
            uKnot = new double[branch.N - uDdim + 1 + uDdim * 2];
            for (int i = 0; i < uDdim; i++)
            {
                uKnot[i] = 0;
            }
            for (int i = 0; i < branch.N - uDdim + 1; i++)
            {
                uKnot[i + uDdim] = i;
            }
            for (int i = 0; i < uDdim; i++)
            {
                uKnot[i + branch.N + 1] = branch.N - uDdim;
            }
            branch.myArch=new Minilla3D.Objects.arch();
            int[] index = new int[uDim];
            for (int i = 1; i < branch.N - uDdim+1; i++)
            {
                for (int l = 0; l < uDim; l++)
                {
                    index[l] = i - 1 + l;
                }
                branch.myArch.elemList.Add(new Minilla3D.Elements.nurbsCurve(uDim, index, i, uKnot));
            }
        }

        void createNurbsElements(leaf leaf)
        {
            double[] uKnot;
            double[] vKnot;

            int N = leaf.nU * leaf.nV;
            int uDim = leaf.srf.OrderU;
            int vDim = leaf.srf.OrderV;
            int uDdim = leaf.srf.OrderU - 1;
            int vDdim = leaf.srf.OrderV - 1;


            uKnot = new double[leaf.nU - uDdim + 1 + uDdim * 2];
            vKnot = new double[leaf.nV - vDdim + 1 + vDdim * 2];
            for (int i = 0; i < uDdim; i++)
            {
                uKnot[i] = 0;
            }
            for (int i = 0; i < vDdim; i++)
            {
                vKnot[i] = 0;
            }
            for (int i = 0; i < leaf.nU - uDdim + 1; i++)
            {
                uKnot[i + uDdim] = i;
            }
            for (int i = 0; i < leaf.nV - vDdim + 1; i++)
            {
                vKnot[i + vDdim] = i;
            }
            for (int i = 0; i < uDdim; i++)
            {
                uKnot[i + leaf.nU + 1] = leaf.nU - uDdim;
            }
            for (int i = 0; i < vDdim; i++)
            {
                vKnot[i + leaf.nV + 1] = leaf.nV - vDdim;
            }
            leaf.myMasonry = new Minilla3D.Objects.masonry();
            for (int j = 1; j < leaf.nV - vDdim + 1; j++)
            {
                for (int i = 1; i < leaf.nU - uDdim + 1; i++)
                {
                    int[] index = new int[uDim * vDim];
                    for (int k = 0; k < vDim; k++)
                    {
                        for (int l = 0; l < uDim; l++)
                        {
                            index[k * uDim + l] = (j - 1 + k) * leaf.nU + i - 1 + l;
                        }
                    }
                    leaf.myMasonry.elemList.Add(new Minilla3D.Elements.nurbsElement(uDim, vDim, index, i, j, uKnot, vKnot));
                }
            }
        }
    }
}
