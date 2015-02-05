using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using ShoNS.Array;
using System.IO;
using System.Reflection;
//using System.Reactive.Linq;
namespace Minilla3D
{
	namespace Objects{
        
		public interface iObject
		{
		}
        public class arch:iObject
        {
            public List<Minilla3D.Elements.nurbsCurve> elemList = new List<Elements.nurbsCurve>();
            public void Add(Minilla3D.Elements.nurbsCurve e)
            {
                elemList.Add(e);
            }
            public void Clear()
            {
                elemList.Clear();
            }
            public void computeGlobalCoord()
            {
                Parallel.ForEach(elemList, (e) =>
                    e.computeGlobalCoord()
                );
            }
            public void setupNodesFromList(double[,] x)
            {
                Parallel.ForEach(elemList, (e) =>
                    e.setupNodesFromList(x)
                    );
            }
            public void setupAiryPotentialFromList(double[] x)
            {
                Parallel.ForEach(elemList, (e) =>
                    e.setupAiryPotentialFromList(x)
                    );
            }
        }
        public class masonry:iObject{
            public List<Minilla3D.Elements.nurbsElement> elemList = new List<Elements.nurbsElement>();
            public void memoryMetric()
            {
                foreach (Minilla3D.Elements.element e in elemList)
                {
                    e.memoryMetric();       //metric->refMetric,invMetric->invRefMetric,dv->refDv
                }
            }
            public void Add(Minilla3D.Elements.nurbsElement e)
            {
                elemList.Add(e);
            }
            public void Clear()
            {
                elemList.Clear();
            }
            public void precompute()
            {
                foreach (var e in elemList)
                {
                    e.precompute();
                }
            }
            public void computeEigenVectors()
            {
                for(int i=0;i<elemList.Count;i++)
                {
                    var e = elemList[i];
                    e.computeEigenVectors();
                }
            }
            public void setupNodesFromList(double[,] x)
            {
                Parallel.ForEach(elemList, (e) =>
                    e.setupNodesFromList(x)
                    );
            }
            public void setupAiryPotentialFromList(double[] x)
            {
                Parallel.ForEach(elemList, (e) =>
                    e.setupAiryPotentialFromList(x)
                    );
            }            
            public void computeHessian()
            {
                Parallel.ForEach(elemList, (e) =>
                    e.computeHessian()
                );
            }
            public void computeGlobalCoord()
            {
                Parallel.ForEach(elemList, (e) =>
                    e.computeGlobalCoord()
                );
            }
            public void computeAngle()
            {
                foreach (var e in elemList)
                {
                    e.computeAngle();
                }
            }
            public int totalNumberOfIntPnts()
            {
                int num = 0;
                foreach (var e in elemList)
                {
                    num += e.nIntPoint;
                }
                return num;
            }
        }

	}
}
