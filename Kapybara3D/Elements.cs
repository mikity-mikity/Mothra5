using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Minilla3D.Elements
{
    public class element
	{
        enum type{
            Cauchy,SPK
        };
        type typeOfStress;
        protected const int __DIM=3;
        public int nIntPoint,nBIntPoint,nNode,elemDim,nDV;                
		public Minilla3D.Elements.integratingPoint[] intP;     //Integrating points
        public Minilla3D.Elements.integratingPoint[] bIntP;     //Integrating Points on border
        protected double[] node;							        //Nodal coordinate (global)
        protected double[] phi;
        protected int[] index;                                    //indeces of the nodes
        double[] gradient;                              //internal force(equivalent nodal force of stress field)
        protected double[,] hess;                                 //Hessian  (Geometric Stiffness only)
        double[] force;                                 //external force(equivalent nodal force of gravity)
        public void giveTension(double t)
        {
            foreach (var p in intP)
            {
                p.giveTension(t);
            }
        }
		public element(int _nNode,int _elemDim,int _nIntPoint)
		{
            var _index=new int[nDV];
			for(int i=0;i<nDV;i++)
			{
				_index[i]=i;
			}
            initialize(_index, _nNode, _elemDim, _nIntPoint);
		}
        public element(int[] _index,int _nNode,int _elemDim,int _nIntPoint)
		{
            initialize(_index,_nNode,_elemDim,_nIntPoint);
		}
        private void initialize(int[] _index,int _nNode,int _elemDim,int _nIntPoint)
        {
            typeOfStress = type.Cauchy;   //default value
            nIntPoint=_nIntPoint;
            nBIntPoint = 0;
            nNode=_nNode;
            elemDim=_elemDim;
            nDV=nNode*__DIM;
            intP=new integratingPoint[nIntPoint];
            this.bIntP = new integratingPoint[0];
            for (int i = 0; i < nIntPoint; i++)
            {
                intP[i]=new integratingPoint(nNode,elemDim);
            }
            node=new double[nDV];
            phi = new double[nNode];
            index=new int[nNode];
		    gradient=new double[nDV];
		    hess=new double [nDV,nDV];
		    force=new double[nDV];
            setupIndex(_index);
        }
        public void setupIndex(int[] _index)
		{
            Array.Copy(_index,index,nNode);
		}
        public void getConnectionCoeff(int i,ref double[,,] _Gamma)
        {
            for (int s = 0; s < 2; s++)
            {
                for (int t = 0; t < 2; t++)
                {
                    for (int u = 0; u < 2; u++)
                    {
                        _Gamma[s, t, u] = this.intP[i].Gamma[s, t, u];
                    }
                }
            }
        }
        public void getBaseVectors(int i, ref double[][] b)
        {
            b[0][0] = this.intP[i].baseVectors[0, 0];
            b[0][1] = this.intP[i].baseVectors[0, 1];
            b[0][2] = this.intP[i].baseVectors[0, 2];
            b[1][0] = this.intP[i].baseVectors[1, 0];
            b[1][1] = this.intP[i].baseVectors[1, 1];
            b[1][2] = this.intP[i].baseVectors[1, 2];
            ///intP[i].f is equivalent to intP[i].baseVectors
            ///intP[i].f is computed in precompute()
            ///intP[i].baseVectos is computed i computeBaseVectors()
        }
        public double[] getNode(int i)
		{
			return new double[3]{node[i*__DIM+0],node[i*__DIM+1],node[i*__DIM+2]};
	    }
		public int getIndex(int i)
		{
			return index[i];
		}
		public void setupNodes(double[,] x)
		{
            Array.Copy(x,node,nDV);
		}
        public void setupAiryPotentialFromList(double[] x)
        {
            for (int i = 0; i < nNode; i++)
            {
                phi[i] = x[index[i]];
            }
        }
        public void setupNodesFromList(double[,] x)
        {
            for (int i = 0; i < nNode; i++)
            {
                for (int j = 0; j < __DIM; j++)
                {
                    node[i * __DIM + j] = x[index[i], j];
                }
            }
        }

		public double Volume{
	        get;
            protected set;
        }
		public double refVolume{
	        get;
            protected set;
        }
        public void computeAiryFunction()
        {
            for (int i = 0; i < nIntPoint; i++)
            {
                intP[i].computeAiryFunction(node);
            }
            for (int i = 0; i < nBIntPoint; i++)
            {
                bIntP[i].computeAiryFunction(node);
            }
        }
        public void precompute()
        {
            for (int i = 0; i < nIntPoint; i++)
            {
                intP[i].precompute(node);
            }
            for (int i = 0; i < nBIntPoint; i++)
            {
                bIntP[i].precompute(node);
            }
        }
        public void computeEigenVectors()
        {
            for(int i=0;i<nIntPoint;i++)
            {
                intP[i].computeEigenVectors();
            }
            for (int i = 0; i < nBIntPoint; i++)
            {
                bIntP[i].computeEigenVectors();
            }
        }
        public void getEigenVectors(double[][] vec, double[] val, int num)
        {
            for (int i = 0; i < 2; i++)
            {
                vec[i][0] = intP[num].eigenVectors[i][0];
                vec[i][1] = intP[num].eigenVectors[i][1];
                vec[i][2] = intP[num].eigenVectors[i][2];
                val[i] = intP[num].eigenValues[i];
            }
        }
        public void getBEigenVectors(double[][] vec, double[] val, int num)
        {
            for (int i = 0; i < 2; i++)
            {
                vec[i][0] = bIntP[num].eigenVectors[i][0];
                vec[i][1] = bIntP[num].eigenVectors[i][1];
                vec[i][2] = bIntP[num].eigenVectors[i][2];
                val[i] = bIntP[num].eigenValues[i];
            }
        }
        public double[] getIntPoint(int i)
        {
            return intP[i].globalCoord;
        }
        public double[] getBIntPoint(int i)
        {
            return bIntP[i].globalCoord;
        }
        public void computeGlobalCoord()
		{
            for (int i = 0; i < nIntPoint; i++)
            {
                intP[i].computeGlobalCoord(node);
            }
            for (int i = 0; i < nBIntPoint; i++)
            {
                bIntP[i].computeGlobalCoord(node);
            }
        }
		public void computeMetric(){
			for(int i=0;i<nIntPoint;i++)
			{
				intP[i].computeMetric(node);
			}
		}
		public void computeBaseVectors()
		{
			for(int i=0;i<nIntPoint;i++)
			{
				intP[i].computeBaseVectors(node);
			}
		}
		public double computeVolume(){
			double v=0;
			for(int i=0;i<nIntPoint;i++)
			{
				v+=intP[i].weight*intP[i].dv;
			}
			this.Volume=v;
			return v;
		}
		public void memoryVolume(){
			this.refVolume=this.Volume;
		}
		public void setRefVolume(double v){
			this.refVolume=v;
		}
		public void memoryMetric(){
			for(int i=0;i<nIntPoint;i++)
			{
				intP[i].memoryMetric();
			}
		}
		public virtual void computeHessian()
		{
			for(int i=0;i<nDV;i++)
			{
				for(int j=0;j<nDV;j++)
				{
					hess[i,j]=0;
				}
			}
			for(int i=0;i<nIntPoint;i++)
			{
				double val=intP[i].weight*intP[i].refDv*0.5;
				for(int j=0;j<elemDim;j++)
				{
				    for(int jj=0;jj<elemDim;jj++)
				    {
				        for(int k=0;k<nDV;k++)
				        {
					        for(int l=0;l<nDV;l++)
					        {
						        double D=intP[i].SPK[j,jj];
						        double S=intP[i].B[j,jj,k,l];
						        hess[k,l]+=val*D*S;
					        }
				        }
                    }
				}
			}
		}
        public void computeEdgeForce(double[,] Force)
        {

            for (int i = 0; i < nDV; i++)
            {
                force[i] = 0;
            }
            for (int i = 0; i < nBIntPoint; i++)
            {
                var p = bIntP[i];
                p.computeEdgeForce();
            }
            //Gauss-Legendre integral
            for (int i = 0; i < nBIntPoint; i++)
            {
                var p = bIntP[i];
                for(int j=0;j<nDV;j++)
                {
                    force[j] += p.weight * p.Fx * p.N[0, j];
                    force[j] += p.weight * p.Fy * p.N[1, j];
                    force[j] += p.weight * p.Fz * p.N[2, j];
                }
            }
            for (int i = 0; i < nNode; i++)
            {
                for (int k = 0; k < 3; k++)
                {
                    Force[index[i], 0] += force[i * 3 + 0];
                    Force[index[i], 1] += force[i * 3 + 1];
                    Force[index[i], 2] += force[i * 3 + 2];
                }
            }
        }
        public void computeGradient()
		{
			
			for(int i=0;i<nDV;i++)
			{
				gradient[i]=0;
			}
            if(typeOfStress==type.Cauchy){
			    for(int i=0;i<nIntPoint;i++)
			    {
				    double val=intP[i].weight*intP[i].dv*0.5;
				    for(int j=0;j<elemDim;j++)
				    {
                        for(int jj=0;jj<elemDim;jj++)
                        {
					        for(int k=0;k<nDV;k++)
					        {
						        for(int l=0;l<nDV;l++)
						        {
							        double D=intP[i].Cauchy[j,jj];
							        double S=intP[i].B[j,jj,k,l];
							        double E=node[l];
							        gradient[k]+=val*D*S*E;
						        }
					        }
                        }
				    }
			    }
            }else{
			    for(int i=0;i<nIntPoint;i++)
			    {
				    double val=intP[i].weight*intP[i].refDv*0.5;
				    for(int j=0;j<elemDim;j++)
				    {
                        for(int jj=0;jj<elemDim;jj++)
                        {
					        for(int k=0;k<nDV;k++)
					        {
						        for(int l=0;l<nDV;l++)
						        {
							        double D=intP[i].SPK[j,jj];
							        double S=intP[i].B[j,jj,k,l];
							        double E=node[l];
							        gradient[k]+=val*D*S*E;
						        }
					        }
                        }
				    }
			    }
            }
		}
		public void mergeHessian(ShoNS.Array.SparseDoubleArray _hess)
		{
			for(int i=0;i<nNode;i++)
			{
				for(int j=0;j<__DIM;j++)
				{
					for(int k=0;k<nNode;k++)
					{
						for(int l=0;l<__DIM;l++)
						{
                            _hess[this.index[i] * __DIM + j, this.index[k] * __DIM + l] += this.hess[i * __DIM + j, k * __DIM + l];
						}
					}
				}
			}
		}
	}
}
