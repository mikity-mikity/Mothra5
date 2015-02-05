using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Minilla3D.Elements
{
    public class nurbsCurve:element
    {
        public class dl
        {
            public nurbsElement.tuple left, right;
            public nurbsElement.tuple target;
            public int index;  //element index
            public int[] internalIndex;
            public int nNode;
            public int elemDim;
            public double area;

            public double ot, t;   //scaled coordinate, coordinate on Rhino, local coordinate on an element
            public double lo;        //local coordinate
            public double[,][] d2;   //gradient of second derivative
            public double[][] d1;    //gradient of first derivative
            public double[] d0;      //gradient of shape function
            public double x, y, z;
            public double[][] gi;
            public double[][] Gi;
            public double[,] gij;
            public double[,] Gij;
            public double[,][] second;
            public double[,] shape;
            public double[, ,] C;
            public double[, , ,] B;
            public double[, , ,] D;
            public double[,] H;
            public double[,] SPK;
            public double dv, refDv;
            public void CtoB(int elemDim, int nDV)
            {
                for (int i = 0; i < elemDim; i++)
                {
                    for (int j = 0; j < elemDim; j++)
                    {
                        for (int u = 0; u < nDV; u++)
                        {
                            for (int v = 0; v < nDV; v++)
                            {
                                B[i, j, u, v] = 0;
                                for (int r = 0; r < __DIM; r++)
                                {
                                    //対称化
                                    B[i, j, u, v] += C[i, r, u] * C[j, r, v];
                                    B[i, j, u, v] += C[j, r, u] * C[i, r, v];
                                }
                            }
                        }
                    }
                }
            }

            public dl(double _ot, double _t, int _index, double _lo, double _area)
            {
                ot = _ot;
                t = _t;
                lo=_lo;

                index = _index;
                area = _area;
                x = 0;
                y = 0;
                z = 0;
                gi = new double[1][] { new double[3] };
                Gi = new double[1][] { new double[3] };
                gij = new double[1, 1];
                Gij = new double[1, 1];
                second = new double[1, 1][] { { new double[3] } };
                H = new double[1, 1];
                SPK = new double[1, 1];
            }
        }
        public void precompute(dl tup)
        {
            //Assume M[0] and M[1] are precomputed,

            //double[][,] M=new double[2][,];
            //M[0]=fM(uNum,_uDim,_uDim-1,uKnot);
            //M[1]=fM(vNum,_vDim,_vDim-1,vKnot);
            tup.internalIndex = this.index;
            tup.nNode = nDV / 3;
            tup.elemDim = elemDim;
            tup.shape = new double[__DIM, nDV];                        //Global coordinate *coefficient*
            tup.C = new double[elemDim, __DIM, nDV];                //Base vectors *coefficient*
            tup.B = new double[elemDim, elemDim, nDV, nDV];          //Metric *coefficient*
            tup.D = new double[elemDim, elemDim, __DIM, nDV];   //Hessian coefficient
            tup.d0 = new double[nDV / 3];
            tup.d1 = new double[elemDim][];
            tup.d2 = new double[elemDim, elemDim][];
            for (int i = 0; i < elemDim; i++)
            {
                tup.d1[i] = new double[nDV / 3];
            }
            for (int i = 0; i < elemDim; i++)
            {
                for (int j = 0; j < elemDim; j++)
                {
                    tup.d2[i, j] = new double[nDV / 3];
                }
            }

            //Shape functions  [N] (for global coordinate)
            double t = tup.lo;
            for (int k = 0; k < dim; k++)
            {
                hh[k] = Math.Pow(t, (dim - k - 1));
            }
            for (int k = 0; k < dim; k++)
            {
                double val = 0;
                for (int l = 0; l < dim; l++)
                {
                    val += hh[l] * M[l, k];
                }
                tt[k] = val;
            }
            for (int j = 0; j < __DIM; j++)
            {
                for (int k = 0; k < nDV; k++)
                {
                    tup.shape[j, k] = 0;
                }
            }
            for (int k = 0; k < nNode; k++)
            {
                //Shape functinos
                double shape = 1.0;
                shape *= tt[dd[k]];
                for (int j = 0; j < __DIM; j++)
                {
                    tup.shape[j, k * __DIM + j] = shape;
                }
            }
            //Create [C]  (for base vectors)
            t = tup.lo;
            {
                for (int k = 0; k < dim - 1; k++)
                {
                    hh[k] = (dim - k - 1) * Math.Pow(t, (dim - k - 2));
                }
                hh[dim - 1] = 0;
            }
            for (int k = 0; k < dim; k++)
            {
                double val = 0;
                for (int l = 0; l < dim; l++)
                {
                    val += hh[l] * M[l, k];
                }
                tt[k] = val;
            }
            for (int jj = 0; jj < __DIM; jj++)
            {
                for (int j = 0; j < nDV; j++)
                {
                    tup.C[0, jj, j] = 0;
                }
            }
            for (int k = 0; k < nNode; k++)
            {
                //[C]
                double C = 1.0;
                for (int j = 0; j < elemDim; j++)
                {
                    C *= tt[dd[k]];
                }
                for (int j = 0; j < __DIM; j++)
                {
                    tup.C[0, j, k * __DIM + j] = C;
                }
            }
            //Create [B]  (for metric)
            tup.CtoB(elemDim, nDV);
            //Create [D] (for second derivative)
            t = tup.lo;
            for (int k = 0; k < dim - 1; k++)
            {
                hh[k] = (dim - k - 1) * (dim - k - 2) * Math.Pow(t, (dim - k - 3));
            }
            hh[dim - 1] = 0;
            hh[dim - 2] = 0;

            for (int k = 0; k < dim; k++)
            {
                double val = 0;
                for (int l = 0; l < dim; l++)
                {
                    val += hh[l] * M[l, k];
                }
                tt[k] = val;
            }
            for (int jj = 0; jj < __DIM; jj++)
            {
                for (int j = 0; j < nDV; j++)
                {
                    tup.D[0, 0, jj, j] = 0;
                }
            }
            for (int k = 0; k < nNode; k++)
            {
                //[D]
                double D = 1.0;
                D *= tt[dd[k]];
                for (int j = 0; j < __DIM; j++)
                {
                    tup.D[0, 0, j, k * __DIM + j] = D;
                }
            }
            compute(tup);
        }
        public void compute(dl tup)
        {
            //Global position
            double X = 0, Y = 0, Z = 0;
            for (int i = 0; i < nDV; i++)
            {
                X += tup.shape[0, i] * node[i];
                Y += tup.shape[1, i] * node[i];
                Z += tup.shape[2, i] * node[i];
            }
            tup.x = X;
            tup.y = Y;
            tup.z = Z;
            //covariant base vectors
            double fx = 0, fy = 0;
            for (int i = 0; i < nDV; i++)
            {
                fx += tup.C[0, 0, i] * node[i];
                fy += tup.C[0, 1, i] * node[i];
            }
            tup.gi[0][0] = fx;
            tup.gi[0][1] = fy;
            tup.gi[0][2] = 0;
            tup.gij[0, 0] = tup.gi[0][0] * tup.gi[0][0] + tup.gi[0][1] * tup.gi[0][1] + tup.gi[0][2] * tup.gi[0][2];
            if (elemDim == 1)
            {
                _inv1(tup.gij, tup.Gij);
            }
            else if (elemDim == 2)
            {
                _inv2(tup.gij, tup.Gij);
            }
            else if (elemDim == 3)
            {
                _inv3(tup.gij, tup.Gij);
            }
            if (elemDim == 1)
            {
                tup.dv = Math.Sqrt(_det1(tup.gij));
            }
            else if (elemDim == 2)
            {
                tup.dv = Math.Sqrt(_det2(tup.gij));
            }
            else if (elemDim == 3)
            {
                tup.dv = Math.Sqrt(_det3(tup.gij));
            }
            tup.refDv = tup.dv;

            //contravatiant base vectors
            double Fx = 0, Fy = 0;
            Fx += tup.gi[0][0] * tup.Gij[0, 0];
            Fy += tup.gi[0][1] * tup.Gij[0, 0];
            tup.Gi[0][0] = Fx;
            tup.Gi[0][1] = Fy;
            tup.Gi[0][2] = 0;
            //Create gradient of hessian with computed connection coefficients
            for (int k = 0; k < nNode; k++)
            {
                tup.d2[0, 0][k] = tup.D[0, 0, 2, k * __DIM + 2];
            }
            for (int k = 0; k < nNode; k++)
            {
                tup.d1[0][k] = tup.C[0, 2, k * __DIM + 2];
            }
            for (int k = 0; k < nNode; k++)
            {
                tup.d0[k] = tup.shape[2, k * __DIM + 2];
            }
        }
        double[,] M;
        int dim;
        double[] hh;
        double[] tt;
        int[] dd;

        double[] _cu;
		double[] _pu;
        public int uDim,vDim;
        private void _inv1(double[,] from, double[,] to)
        {
            if (from[0, 0] <= 0)
            {
                to[0, 0] = 0;
                from[0, 0] = 0;
            }
            else
            {
                to[0, 0] = 1 / from[0, 0];
            }
        }
        private void _inv2(double[,] from, double[,] to)
        {
            double det = _det2(from);
            if (det <= 0)
            {
                to[0, 0] = 0;
                to[1, 0] = 0;
                to[0, 1] = 0;
                to[1, 1] = 0;
                from[0, 0] = 0;
                from[1, 0] = 0;
                from[0, 1] = 0;
                from[1, 1] = 0;
            }
            else
            {
                to[0, 0] = from[1, 1] / det;
                to[1, 1] = from[0, 0] / det;
                to[1, 0] = -from[1, 0] / det;
                to[0, 1] = -from[0, 1] / det;
            }
        }
        private void _inv3(double[,] from, double[,] to)
        {
            double det = _det3(from);
            if (det <= 0)
            {
                to[0, 0] = 0;
                to[0, 1] = 0;
                to[0, 2] = 0;
                to[1, 0] = 0;
                to[1, 1] = 0;
                to[1, 2] = 0;
                to[1, 3] = 0;
                to[2, 0] = 0;
                to[2, 1] = 0;
                to[2, 2] = 0;
                from[0, 0] = 0;
                from[0, 1] = 0;
                from[0, 2] = 0;
                from[1, 0] = 0;
                from[1, 1] = 0;
                from[1, 2] = 0;
                from[1, 3] = 0;
                from[2, 0] = 0;
                from[2, 1] = 0;
                from[2, 2] = 0;
            }
            else
            {
                to[0, 0] = (from[1, 1] * from[2, 2] - from[1, 2] * from[2, 1]) / det;
                to[1, 0] = (from[2, 1] * from[0, 2] - from[2, 2] * from[0, 1]) / det;
                to[2, 0] = (from[0, 1] * from[1, 2] - from[0, 2] * from[1, 1]) / det;
                to[0, 1] = (from[1, 2] * from[2, 0] - from[1, 0] * from[2, 2]) / det;
                to[1, 1] = (from[2, 2] * from[0, 0] - from[2, 0] * from[0, 2]) / det;
                to[2, 1] = (from[0, 2] * from[1, 0] - from[0, 0] * from[1, 2]) / det;
                to[0, 2] = (from[1, 0] * from[2, 1] - from[1, 1] * from[2, 0]) / det;
                to[1, 2] = (from[2, 0] * from[0, 1] - from[2, 1] * from[0, 0]) / det;
                to[2, 2] = (from[0, 0] * from[1, 1] - from[0, 1] * from[1, 0]) / det;
            }
        }
        private double _det1(double[,] m)
        {
            return m[0, 0];
        }
        private double _det2(double[,] m)
        {
            return m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0];
        }
        private double _det3(double[,] m)
        {
            return m[0, 0] * m[1, 1] * m[2, 2]
            + m[1, 0] * m[2, 1] * m[0, 2]
            + m[2, 0] * m[0, 1] * m[1, 2]
            - m[2, 0] * m[1, 1] * m[0, 2]
            - m[1, 0] * m[0, 1] * m[2, 2]
            - m[0, 0] * m[2, 1] * m[1, 2];
        }
        double[] fN(int _i, int _k, int _dim, int dim, double[] knot)
		{
		    if (_dim==1)
			{
		        double[] F=new double[dim]; 
				for (int i=0;i<dim;i++)
				{
					F[i]=0;
				}
		        if (_k==_i)
				{
					F[dim-1]=1;
				}
		        return F;
			}
		    double[] S1=fN(_i,_k,_dim-1,dim,knot);
			double[] S2=fN(_i,_k+1,_dim-1,dim,knot);
			double E1=knot[_k+_dim-2]-knot[_k-1];
			double E2=knot[_k+_dim-1]-knot[_k];
			double[] D1=new double[2]{0,0};
			double[] D2=new double[2]{0,0};
		    if (E1>0)
			{
				D1[0]=1d/E1;
				D1[1]=-knot[_k-1]/E1;
			}
			if (E2>0)
			{
				D2[0]=-1d/E2;
				D2[1]=knot[_k+_dim-1]/E2;
			}
		    double[] F2=new double[dim]; 
			for (int i=0;i<dim;i++)
			{
				F2[i]=0;
			}
			for(int i=1;i<dim;i++)
			{
				F2[i-1]=F2[i-1]+S1[i]*D1[0];
				F2[i]=F2[i]+S1[i]*D1[1];
				F2[i-1]=F2[i-1]+S2[i]*D2[0];
				F2[i]=F2[i]+S2[i]*D2[1];
			}
			return F2;
		}
        double[,] fM(int shift, int dim, int ddim, double[] knot){
			double[,] M=new double[dim,dim];
            for (int i = 0; i < dim; i++)
            {
                for (int j = 0; j < dim; j++)
                {
                    M[i,j] = 0;
                }
            }
		    for(int k=shift;k<dim+shift;k++)
			{
		        double[] D=fN(shift+ddim,k,dim,dim,knot);
				for (int n =0;n<dim;n++)
				{
					M[n,k-shift]=D[n];
				}
			}

			double[,] S=new double[dim,dim];
            for(int i=0;i<dim;i++)
			{
                for(int j=0;j<dim;j++)
			    {
				    S[i,j]=0;
			    }
            }
			for(int  n =1;n<dim+1;n++)
			{
				for (int k=1+n;k<dim+2;k++)
				{
					if (n==dim)
					{
						for (int t =0;t<n-1;t++)
						{
						   S[t,n-1]=0;
						}
						S[(n-1),n-1]=1;
					}else
					{
						S[(k-2),n-1]=binominal(dim-n,dim+1-k)*Math.Pow(shift-1,k-1-n);
					}
				}
			}
			double[,] G=new double[dim,dim];
			for (int j=0;j<dim;j++)
			{
				for(int k=0;k<dim;k++)
				{
					double v=0;
					for(int l=0;l<dim;l++)
					{
						v+=S[j,l]*M[l,k];
					}
					G[j,k]=v;
				}
			}
			return G;
		}
        private static int Factorial(int x)
        {
            if (x == 0) return 1;
            if (x == 1) return 1;
            if (x == 2) return 2;
            if (x == 3) return 6;
            if (x == 4) return 24;
            if (x == 5) return 120;
            int val = 1;
            for (int i = 2; i <= x; i++)
            {
                val *= i;
            }
            return val;
        }

        public static double binominal(int N, int k)
        {
            return Factorial(N) / Factorial(N - k) / Factorial(k);
        }

        public void setPlane(double a, double b, double c, double d)
        {
            foreach (var p in intP)
            {
                p.refIntP.setPlane(a, b, c, d);
            }
        }
        public int numberOfConstraintConditions()
        {
            return nIntPoint;
        }
        public int mergeResidual(ShoNS.Array.DoubleArray residual, int i)
        {
            for (int j = 0; j < nBIntPoint; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    double resid = bIntP[j].getResidualOfBoundaryCondition(node,k);
                    residual[i + j * 2 + k] = resid;
                }
            }
            return i + nBIntPoint * 2;
        }
        public int mergeJacobian(ShoNS.Array.SparseDoubleArray jacobian, int i)
        {
            for (int j = 0; j < nBIntPoint; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    var grad = bIntP[j].getGradientOfBoundaryCondition(k);
                    for(int f=0;f<nNode;f++)
                    {
                        jacobian[i + j*2 + k, index[f]] = grad[f];
                    }
                }
            }
            return i + nBIntPoint * 2;
        }
        public int mergeJacobianOfPosition(ShoNS.Array.SparseDoubleArray jacob, int num)
        {
            for (int i = 0; i < nIntPoint; i++)
            {
                var grad = intP[i].N;
                for (int j = 0; j < nNode; j++)
                {
                    jacob[num + i, index[j]] = grad[2, j * 3 + 2];
                }
            }
            return num + nIntPoint;
        }

        double __cu(int n, int dim)
        {
            switch (dim)
            {
                case 6:
                    switch (n)
                    {
                        case 0:
                            return (-0.9491079123427585245261897) * 0.5 + 0.5;
                        case 1:
                            return (-0.7415311855993944398638648) * 0.5 + 0.5;
                        case 2:
                            return (-0.4058451513773971669066064) * 0.5 + 0.5;
                        case 3:
                            return 0.5;
                        case 4:
                            return (0.4058451513773971669066064) * 0.5 + 0.5;
                        case 5:
                            return (0.7415311855993944398638648) * 0.5 + 0.5;
                        case 6:
                            return (0.9491079123427585245261897) * 0.5 + 0.5;
                        default:
                            return 0;
                    }
                case 5:
                    switch (n)
                    {
                        case 0:
                            return (-0.9324695142031521) * 0.5 + 0.5;
                        case 1:
                            return (-0.6612093864662645) * 0.5 + 0.5;
                        case 2:
                            return (-0.2386191860831969) * 0.5 + 0.5;
                        case 3:
                            return (0.2386191860831969) * 0.5 + 0.5;
                        case 4:
                            return (0.6612093864662645) * 0.5 + 0.5;
                        case 5:
                            return (0.9324695142031521) * 0.5 + 0.5;
                        default:
                            return 0;
                    }
                case 4:
                    switch (n)
                    {
                        case 0:
                            return (-0.9061798459386640)*0.5+0.5;
                        case 1:
                            return (-0.5384693101056831)*0.5+0.5;
                        case 2:
                            return (0.0000000000000000)*0.5+0.5;
                        case 3:
                            return (0.5384693101056831)*0.5+0.5;
                        case 4:
                            return (0.9061798459386640)*0.5+0.5;
				        default:
                            return 0;
                    }
                case 3:
                    switch (n)
                    {
                        case 0:
                            return (-0.8611363115940526)*0.5+0.5;
                        case 1:
                            return (-0.3399810435848563)*0.5+0.5;
                        case 2:
                            return (0.3399810435848563)*0.5+0.5;
                        case 3:
                            return (0.8611363115940526)*0.5+0.5;
                        default:
                            return 0;
                    }
                case 2:
                    switch (n)
                    {
                        case 0:
                            return (-0.7745966692414834) * 0.5 + 0.5;
                        case 1:
                            return (0.0000000000000000) * 0.5 + 0.5;
                        case 2:
                            return (0.7745966692414834) * 0.5 + 0.5;
                        default:
                            return 0;
                    }
                default:
                    return 0;

            }
        }
        double __pu(int n, int dim)
        {
            switch (dim)
            {
                case 6:
                    switch (n)
                    {
                        case 0:
                            return 0.1294849661688696932706114 * 0.5;
                        case 1:
                            return 0.2797053914892766679014678 * 0.5;;
                        case 2:
                            return 0.3818300505051189449503698 * 0.5;;
                        case 3:
                            return 0.4179591836734693877551020 * 0.5;;
                        case 4:
                            return 0.3818300505051189449503698 * 0.5;;
                        case 5:
                            return 0.2797053914892766679014678 * 0.5;;
                        case 6:
                            return 0.1294849661688696932706114 * 0.5;;
                        default:
                            return 0;
                    }
                case 5:
                    switch (n)
                    {
                        case 0:
                            return 0.1713244923791704 * 0.5;
                        case 1:
                            return 0.3607615730481386 * 0.5;
                        case 2:
                            return 0.4679139345726910 * 0.5;
                        case 3:
                            return 0.4679139345726910 * 0.5;
                        case 4:
                            return 0.3607615730481386 * 0.5;
                        case 5:
                            return 0.1713244923791704 * 0.5;
                        default:
                            return 0;
                    }
                case 4:
                    switch (n)
                    {
                        case 0:
                            return 0.2369268850561891 * 0.5;
                        case 1:
                            return 0.4786286704993665 * 0.5;
                        case 2:
                            return 0.5688888888888889 * 0.5;
                        case 3:
                            return 0.4786286704993665 * 0.5;
                        case 4:
                            return 0.2369268850561891 * 0.5;
                        default:
                            return 0;
                    }
                case 3:
                    switch (n)
                    {
                        case 0:
                            return 0.3478548451374538*0.5;
                        case 1:
                            return 0.6521451548625461*0.5;
                        case 2:
                            return 0.6521451548625461*0.5;
                        case 3:
                            return 0.3478548451374538 * 0.5;
                        default:
                            return 0;
                    }
                case 2:
                    switch (n)
                    {
                        case 0:
                            return 0.5555555555555556*0.5;
                        case 1:
                            return 0.8888888888888888*0.5;
                        case 2:
                            return 0.5555555555555556 * 0.5;
                        default:
                            return 0;
                    }
                default:
                    return 0;
            }
        }
        public nurbsCurve(int _uDim, int[] _index, int uNum, double[] uKnot)
            : base(_index, _uDim, 1, (_uDim + 1))
        {
            uDim = _uDim;
            _cu = new double[uDim + 1];
            _pu = new double[uDim + 1];
            int[] ss = new int[nIntPoint];		//Indeces for integrating points
			dd=new int[nNode];		    //Indeces for nodes
            dim=uDim;
			hh=new double[uDim];
		    tt=new double[uDim];

					
			//Weight coefficient distribution
			//Coordinates distribution
            for (int i = 0; i < uDim+1; i++)
            {
                _cu[i] = __cu(i, uDim);
                _pu[i] = __pu(i, uDim);
            }
            

            //Indeces for inthegrating points
            ss[0] = 0;
            for (int i = 1; i < nIntPoint; i++)
			{
				ss[i]=ss[i-1];
				if(ss[i]<dim)
				{
					ss[i]++;
				}
			}

			//Indices for nodes
		    dd[0]=0;
			for(int i=1;i<nNode;i++)
			{
			    dd[i]=dd[i-1];
				if(dd[i]<dim-1)
				{
					dd[i]++;
				}
			}
			//weight coefficients and loacl coordinates for integrating points
			for(int i=0;i<nIntPoint;i++)
			{
				intP[i].weight=1.0;
				intP[i].localCoord[0]=_cu[ss[i]];
				intP[i].weight*=_pu[ss[i]];
			}
		    M=fM(uNum,_uDim,_uDim-1,uKnot);
			//Shape functions  [N] (for global coordinate)
			for(int i=0;i<nIntPoint;i++)
			{
				double t=intP[i].localCoord[0];
				for(int k=0;k<dim;k++)
				{
					hh[k]=Math.Pow(t,(dim-k-1));
				}
				for(int k=0;k<dim;k++)
				{
					double val=0;
					for(int l=0;l<dim;l++)
					{
						val+=hh[l]*M[l,k];
					}
					tt[k]=val;
				}
				for(int j=0;j<__DIM;j++)
				{
                    for(int k=0;k<nDV;k++)
                    {
					    intP[i].N[j,k]=0;
                    }
				}
				for(int k=0;k<nNode;k++)
				{
					//Shape functinos
					double N=1.0;
					N*=tt[dd[k]];
					for(int j=0;j<__DIM;j++)
					{
						intP[i].N[j,k*__DIM+j]=N;
					}
				}

                //Create [C]  (for base vectors)
				t=intP[i].localCoord[0];
				{
					for(int k=0;k<dim-1;k++)
					{
						hh[k]=(dim-k-1)*Math.Pow(t,(dim-k-2));
					}
					hh[dim-1]=0;
				}
				for(int k=0;k<dim;k++)
				{
					double val=0;
					for(int l=0;l<dim;l++)
					{
						val+=hh[l]*M[l,k];
					}
					tt[k]=val;
				}
                for(int jj=0;jj<__DIM;jj++)
                {
					for(int j=0;j<nDV;j++)
					{
						intP[i].C[0,jj,j]=0;
					}
                }
				for(int k=0;k<nNode;k++)
				{
					//[C]
					double C=1.0;
					C*=tt[dd[k]];
					for(int j=0;j<__DIM;j++)
					{
						intP[i].C[0,j,k*__DIM+j]=C;
					}
				}
				//Create [B]  (for metric)
                intP[i].CtoB();

                //Create [D] (for second derivative)
                t = intP[i].localCoord[0];
                for (int k = 0; k < dim - 1; k++)
                {
                    hh[k] = (dim - k - 1) *(dim - k - 2) * Math.Pow(t, (dim - k - 3));
                }
                hh[dim - 1] = 0;
                hh[dim - 2] = 0;

                for (int k = 0; k < dim; k++)
                {
                    double val = 0;
                    for (int l = 0; l < dim; l++)
                    {
                        val += hh[l] * M[l, k];
                    }
                    tt[k] = val;
                }
                        
                for (int jj = 0; jj < __DIM; jj++)
                {
                    for (int j = 0; j < nDV; j++)
                    {
                        intP[i].D[0, 0,jj, j] = 0;
                    }
                }
                for (int k = 0; k < nNode; k++)
                {
                    //[D]
                    double D = 1.0;
                    D *= tt[dd[k]];
                    for (int j = 0; j < __DIM; j++)
                    {
                        intP[i].D[0, 0, j, k * __DIM + j] = D;
                    }
                }
            }
        }

    }
}
