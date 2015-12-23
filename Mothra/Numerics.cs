using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ShoNS.Array;
using Rhino.Geometry;
namespace mikity.ghComponents
{
    class msgclass : mosek.Stream
    {
        string prefix;
        public msgclass(string prfx)
        {
            prefix = prfx;
        }

        public override void streamCB(string msg)
        {
            Console.Write("{0}{1}", prefix, msg);
        }
    }

    public partial class Mothra4 : Grasshopper.Kernel.GH_Component
    {
        public void defineKinkAngle(branch branch,leaf leaf,mosek.Task task, int numCon, int numVar)
        {
            //todo add valDc//
            for (int t= 0; t < branch.tuples.Count(); t++)
            {
                task.putaij(numCon + t, numVar + t, -1);
                task.putconbound(numCon + t, mosek.boundkey.fx, 0, 0);
                var tup = branch.tuples[t];
                var target = tup.target;
                target.dcdtstar[0] = target.dcdt[1];
                target.dcdtstar[1] = -target.dcdt[0];
                double gamma = 0;
                for (int i = 0; i < 2; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        gamma += target.dcdt[i] * target.gij[i, j] * target.dcdt[j];
                    }
                }
                for (int k = 0; k < target.nNode; k++)
                {
                    for (int i = 0; i < 2; i++)
                    {
                        target.s[i] =target.d1[i][k];
                    }
                    var val = 0d;
                    for (int i = 0; i < 2; i++)
                    {
                        for (int j = 0; j < 2; j++)
                        {
                            val += target.s[i] * target.Gij[i, j] * target.dcdtstar[j];
                        }
                    }
                    val *= target.refDv;
                    val /= Math.Sqrt(gamma);
                    task.putaij(numCon + t, target.internalIndex[k] + leaf.varOffset, val);
                }
            }
        }
        public void defineKinkAngleC(branch branch, leaf leaf, mosek.Task task, int numCon, int numVar)
        {
            //todo add valDc//
            for (int t = 0; t < branch.tuples.Count(); t++)
            {
                task.putaij(numCon + t, numVar + t, -1);
                task.putconbound(numCon + t, mosek.boundkey.fx, 0, 0);
                var tup = branch.tuples[t];
                var target = tup.target;
                target.dcdtstar[0] = target.dcdt[1];
                target.dcdtstar[1] = -target.dcdt[0];
                double gamma = 0;
                for (int i = 0; i < 2; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        gamma += target.dcdt[i] * target.gij[i, j] * target.dcdt[j];
                    }
                }
                for (int k = 0; k < target.nNode; k++)
                {
                    for (int i = 0; i < 2; i++)
                    {
                        target.s[i] = target.d1[i][k];
                    }
                    var val = 0d;
                    for (int i = 0; i < 2; i++)
                    {
                        for (int j = 0; j < 2; j++)
                        {
                            val += target.s[i] * target.Gij[i, j] * target.dcdtstar[j];
                        }
                    }
                    val *= target.refDv;
                    val /= Math.Sqrt(gamma);
                    task.putaij(numCon + t, target.internalIndex[k] + leaf.varOffset, val);
                }
                //a
                for (int i = 0; i < 2; i++)
                {
                    target.s[i] = target.gi[i][0];
                }
                double val2 = 0d;
                for (int i = 0; i < 2; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        val2 += target.s[i] * target.Gij[i, j] * target.dcdtstar[j];
                    }
                }
                val2 *= target.refDv;
                val2 /= Math.Sqrt(gamma);
                task.putaij(numCon + t, branch.slice.varOffset, val2);

                //b
                for (int i = 0; i < 2; i++)
                {
                    target.s[i] = target.gi[i][1];
                }
                val2 = 0d;
                for (int i = 0; i < 2; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        val2 += target.s[i] * target.Gij[i, j] * target.dcdtstar[j];
                    }
                }
                val2 *= target.refDv;
                val2 /= Math.Sqrt(gamma);
                task.putaij(numCon + t, branch.slice.varOffset + 1, val2);
            }
        }
        public void defineKinkAngle2(branch branch, leaf leaf1, leaf leaf2, mosek.Task task, int numCon, int numVar)
        {
            for (int t = 0; t < branch.tuples.Count(); t++)
            {
                task.putaij(numCon + t, numVar + t, -1);
                task.putconbound(numCon + t, mosek.boundkey.fx, 0, 0);
                var tup = branch.tuples[t];
                for (int h = 0; h < 2; h++)
                {
                    Minilla3D.Elements.nurbsElement.tuple target = null;
                    leaf leaf = null;
                    if (h == 0) { target = tup.left; leaf = leaf1; }
                    if (h == 1) { target = tup.right; leaf = leaf2; }
                    target.dcdtstar[0] = target.dcdt[1];
                    target.dcdtstar[1] = -target.dcdt[0];
                    double gamma = 0;
                    for (int i = 0; i < 2; i++)
                    {
                        for (int j = 0; j < 2; j++)
                        {
                            gamma += target.dcdt[i] * target.gij[i, j] * target.dcdt[j];
                            gamma += target.dcdt[i] * target.gij[i, j] * target.dcdt[j];
                        }
                    }
                    for (int k = 0; k < target.nNode; k++)
                    {
                        for (int i = 0; i < 2; i++)
                        {
                            target.s[i] = target.d1[i][k];
                        }
                        var val = 0d;
                        for (int i = 0; i < 2; i++)
                        {
                            for (int j = 0; j < 2; j++)
                            {
                                val += target.s[i] * target.Gij[i, j] * target.dcdtstar[j];
                            }
                        }
                        val *= target.refDv;
                        val /= Math.Sqrt(gamma);
                        task.putaij(numCon + t, target.internalIndex[k] + leaf.varOffset, val);
                    }
                }
            }
        }

        public void tieBranchD1(branch branch,leaf leaf,mosek.Task task,int num0/*1 or 2*/,int num1/*0 or 1*/)
        {
            if (leaf.branch[0] == branch)
            {
                if (leaf.nU == branch.N)
                {
                    for (int i = 0; i < branch.N; i++)
                    {
                        task.putconbound(i * num0 + num1 + branch.conOffset , mosek.boundkey.fx, 0, 0);
                        if (leaf.flip[0])
                        {
                            task.putaij(branch.conOffset + i * num0 + num1 , i + branch.varOffset, 1);
                            task.putaij(branch.conOffset + i * num0 + num1 , (branch.N - 1 - i) + leaf.varOffset, -1);
                        }
                        else
                        {
                            task.putaij(branch.conOffset + i * num0 + num1 , i + branch.varOffset, 1);
                            task.putaij(branch.conOffset + i * num0 + num1 , i + leaf.varOffset, -1);
                        }
                    }
                }
            }
            if (leaf.branch[1] == branch)
            {
                if (leaf.nV == branch.N)
                {
                    for (int i = 0; i < branch.N; i++)
                    {
                        task.putconbound(i * num0 + num1 + branch.conOffset, mosek.boundkey.fx, 0, 0);
                        if (leaf.flip[1])
                        {
                            task.putaij(branch.conOffset + i * num0 + num1, i + branch.varOffset, 1);
                            task.putaij(branch.conOffset + i * num0 + num1, leaf.nU * (branch.N - i) - 1 + leaf.varOffset, -1);
                        }
                        else
                        {
                            task.putaij(branch.conOffset + i * num0 + num1, i + branch.varOffset, 1);
                            task.putaij(branch.conOffset + i * num0 + num1, leaf.nU * (i + 1) - 1 + leaf.varOffset, -1);
                        }
                    }
                }
            }
            if (leaf.branch[2] == branch)
            {
                if (leaf.nU == branch.N)
                {
                    for (int i = 0; i < branch.N; i++)
                    {
                        task.putconbound(i * num0 + num1 + branch.conOffset, mosek.boundkey.fx, 0, 0);
                        if (leaf.flip[2])
                        {
                            task.putaij(branch.conOffset + i * num0 + num1, i + branch.varOffset, 1);
                            task.putaij(branch.conOffset + i * num0 + num1, leaf.nU * (leaf.nV - 1) + (branch.N - 1 - i) + leaf.varOffset, -1);
                        }
                        else
                        {
                            task.putaij(branch.conOffset + i * num0 + num1, i + branch.varOffset, 1);
                            task.putaij(branch.conOffset + i * num0 + num1, leaf.nU * (leaf.nV - 1) + i + leaf.varOffset, -1);
                        }
                    }
                }
            }

            if (leaf.branch[3] == branch)
            {
                if (leaf.nV == branch.N)
                {
                    for (int i = 0; i < branch.N; i++)
                    {
                        task.putconbound(i * num0 + num1 + branch.conOffset, mosek.boundkey.fx, 0, 0);
                        if (leaf.flip[3])
                        {
                            task.putaij(branch.conOffset + i * num0 + num1, i + branch.varOffset, 1);
                            task.putaij(branch.conOffset + i * num0 + num1, leaf.nU * (branch.N - 1 - i) + leaf.varOffset, -1);
                        }
                        else
                        {
                            task.putaij(branch.conOffset + i * num0 + num1, i + branch.varOffset, 1);
                            task.putaij(branch.conOffset + i * num0 + num1, leaf.nU * i + leaf.varOffset, -1);
                        }
                    }
                }
            }

        }
        public void tieBranchD3(branch branch, leaf leaf, mosek.Task task, int num0/*1 or 2*/, int num1/*0 or 1*/)
        {
            if (leaf.branch[0] == branch)
            {
                if (leaf.nU == branch.N)
                {
                    for (int i = 0; i < branch.N; i++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            task.putconbound((i * num0 + num1) * 3 + k + branch.conOffset, mosek.boundkey.fx, 0, 0);
                            if (leaf.flip[0])
                            {
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (i) * 3 + k + branch.varOffset, 1);
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, ((branch.N - 1 - i) * 3 + k) + leaf.varOffset, -1);
                            }
                            else
                            {
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (i) * 3 + k + branch.varOffset, 1);
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (i) * 3 + k + leaf.varOffset, -1);
                            }
                        }
                    }
                }
            }
            if (leaf.branch[1] == branch)
            {
                if (leaf.nV == branch.N)
                {
                    for (int i = 0; i < branch.N; i++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            task.putconbound((i * num0 + num1) * 3 + k + branch.conOffset, mosek.boundkey.fx, 0, 0);
                            if (leaf.flip[1])
                            {
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (i) * 3 + k + branch.varOffset, 1);
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (leaf.nU * (branch.N - i) - 1) * 3 + k + leaf.varOffset, -1);
                            }
                            else
                            {
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (i) * 3 + k + branch.varOffset, 1);
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (leaf.nU * (i + 1) - 1) * 3 + k + leaf.varOffset, -1);
                            }
                        }
                    }
                }
            }
            if (leaf.branch[2] == branch)
            {
                if (leaf.nU == branch.N)
                {
                    for (int i = 0; i < branch.N; i++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            task.putconbound((i * num0 + num1) * 3 + k + branch.conOffset, mosek.boundkey.fx, 0, 0);
                            if (leaf.flip[2])
                            {
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (i) * 3 + k + branch.varOffset, 1);
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (leaf.nU * (leaf.nV - 1) + (branch.N - 1 - i)) * 3 + k + leaf.varOffset, -1);
                            }
                            else
                            {
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (i) * 3 + k + branch.varOffset, 1);
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (leaf.nU * (leaf.nV - 1) + i) * 3 + k + leaf.varOffset, -1);
                            }
                        }
                    }
                }
            }

            if (leaf.branch[3] == branch)
            {
                if (leaf.nV == branch.N)
                {
                    for (int i = 0; i < branch.N; i++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            task.putconbound((i * num0 + num1)*3+k + branch.conOffset, mosek.boundkey.fx, 0, 0);
                            if (leaf.flip[3])
                            {
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (i) * 3 + k + branch.varOffset, 1);
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (leaf.nU * (branch.N - 1 - i)) * 3 + k + leaf.varOffset, -1);
                            }
                            else
                            {
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (i) * 3 + k + branch.varOffset, 1);
                                task.putaij(branch.conOffset + (i * num0 + num1) * 3 + k, (leaf.nU * i) * 3 + k + leaf.varOffset, -1);
                            }
                        }
                    }
                }
            }

        }
        public void mosek2(List<leaf> _listLeaf, List<branch> _listBranch, List<node> _listNode, double force,bool fix,double exp)
        {
            //variable settings
            ShoNS.Array.SparseDoubleArray mat = new SparseDoubleArray(_listNode.Count * 3, _listNode.Count * 3);
            ShoNS.Array.SparseDoubleArray F = new SparseDoubleArray(_listNode.Count * 3, 1);
            ShoNS.Array.SparseDoubleArray xx = new SparseDoubleArray(_listNode.Count*3,1);
            ShoNS.Array.SparseDoubleArray shift = new SparseDoubleArray(_listNode.Count*3, _listNode.Count*3);
            for (int k = 0; k < _listNode.Count;k++ )
            {
                var node = _listNode[k];
                xx[k * 3 + 0, 0] = node.x;
                xx[k * 3 + 1, 0] = node.y;
                xx[k * 3 + 2, 0] = node.z;
                if (node.nodeType == node.type.fx)
                {
                    for (int i = 0; i < node.shareB.Count; i++)
                    {
                        var branch = node.shareB[i];
                        var index = node.numberB[i];
                        if (branch.branchType == branch.type.fix)
                        {
                            node.x = branch.crv.Points[index].Location.X;
                            node.y = branch.crv.Points[index].Location.Y;
                            node.z = branch.slice2.height;
                            xx[k * 3 + 0, 0] = node.x;
                            xx[k * 3 + 1, 0] = node.y;
                            xx[k * 3 + 2, 0] = node.z;
                        }
                        if (branch.endTypeL == branch.endtype.fix && node.numberB[i] == 0)
                        {
                            node.x = branch.crv.Points[index].Location.X;
                            node.y = branch.crv.Points[index].Location.Y;
                            node.z = branch.sliceL.height;
                            xx[k * 3 + 0, 0] = node.x;
                            xx[k * 3 + 1, 0] = node.y;
                            xx[k * 3 + 2, 0] = node.z;

                        }
                        if (branch.endTypeL == branch.endtype.fix && node.numberB[i] == branch.crv.Points.Count-1)
                        {
                            node.x = branch.crv.Points[index].Location.X;
                            node.y = branch.crv.Points[index].Location.Y;
                            node.z = branch.sliceR.height;
                            xx[k * 3 + 0, 0] = node.x;
                            xx[k * 3 + 1, 0] = node.y;
                            xx[k * 3 + 2, 0] = node.z;

                        }
                    }
                }
            }
            
            List<int> series=new List<int>();
            List<Point3d> origin = new List<Point3d>();
            for(int i=0;i<_listNode.Count;i++)
            {
                series.Add(i);
            }
            int L1=0;
            int L2=_listNode.Count;
            for(int i=0;i<_listNode.Count;i++)
            {
                var node=_listNode[i];
                if(node.nodeType==node.type.fx)
                {
                    L2--;
                    series[i]=L2;
                }else{
                    series[i]=L1;
                    origin.Add(new Point3d(node.x, node.y, node.z));
                    L1++;
                }
            }
            for (int i = 0; i < _listNode.Count; i++)
            {
                shift[i * 3 + 0, series[i] * 3 + 0] = 1;
                shift[i * 3 + 1, series[i] * 3 + 1] = 1;
                shift[i * 3 + 2, series[i] * 3 + 2] = 1;
                //F[i * 3 + 2, 0] = -force;//force
            }
            foreach (var leaf in _listLeaf)
            {
                foreach (var tup in leaf.tuples)
                {
                    var d0=tup.d0;
                    for(int i = 0; i < tup.nNode; i++)
                    {
                        int k = 2;
                        var val = d0[i] * tup.refDv * tup.area*force;
                        //GGSystem.Windows.Forms.MessageBox.Show(tup.eigenValues[0].ToString() + "," + tup.eigenValues[1].ToString());
                        double factor = Math.Abs(Math.Max(tup.eigenValues[0], tup.eigenValues[1]));
                        F[leaf.globalIndex[tup.internalIndex[i]] * 3 + k,0] -= val * Math.Pow(factor, exp);
                    }
                }
            }
            foreach (var leaf in _listLeaf)
            {
                foreach (var tup in leaf.tuples)
                {
                    var det = tup.SPK[0, 0] * tup.SPK[1, 1] - tup.SPK[0, 1] * tup.SPK[0, 1];
                    if (det > 0)
                    {
                        for (int i = 0; i < tup.nNode; i++)
                        {
                            for (int j = 0; j < tup.nNode; j++)
                            {
                                for (int k = 0; k < 3; k++)
                                {
                                    for (int l = 0; l < 2; l++)
                                    {
                                        for (int m = 0; m < 2; m++)
                                        {
                                            var val = tup.B[l, m, i * 3 + k, j * 3 + k] * tup.SPK[l, m] * tup.refDv * tup.area; //tup.area is just a weight of Gauss integration
                                            mat[leaf.globalIndex[tup.internalIndex[i]] * 3 + k, leaf.globalIndex[tup.internalIndex[j]] * 3 + k] += val;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            foreach (var branch in _listBranch)
            {
                foreach (var tup in branch.tuples)
                {
                    for (int i = 0; i < tup.nNode; i++)
                    {
                        for (int j = 0; j < tup.nNode; j++)
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                //if (tup.SPK[0, 0] > 0)
                                {
                                    var val = tup.B[0, 0, i * 3 + k, j * 3 + k] * tup.SPK[0, 0] * tup.refDv * tup.area;
                                    for (int l = 0; l < 1; l++)
                                    {
                                        for (int m = 0; m < 1; m++)
                                        {
                                            var val2 = tup.B[l, m, i * 3 + k, j * 3 + k] * tup.SPK[l, m] * tup.refDv * tup.area;
                                            mat[branch.globalIndex[tup.internalIndex[i]] * 3 + k, branch.globalIndex[tup.internalIndex[j]] * 3 + k] += val2;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            var reactionForce = mat * xx as SparseDoubleArray;
            double max1 = 0;
            double max2 = 0;
            for (int i = 0; i < _listNode.Count; i++)
            {
                if (_listNode[i].nodeType == node.type.fr)
                {
                    var nx = reactionForce[i * 3 + 0, 0];
                    var ny = reactionForce[i * 3 + 1, 0];
                    var norm = Math.Sqrt(nx * nx + ny * ny);
                    if (norm > max1) max1 = norm;
                }
                else
                {
                    var nx = reactionForce[i * 3 + 0, 0];
                    var ny = reactionForce[i * 3 + 1, 0];
                    var norm = Math.Sqrt(nx * nx + ny * ny);
                    if (norm > max2) max2 = norm;
                }
            }
            //System.Windows.Forms.MessageBox.Show(max1.ToString() + ".-." + max2.ToString());
            var newMat = (shift.T.Multiply(mat) as SparseDoubleArray).Multiply(shift) as SparseDoubleArray;
            var newxx = shift.T.Multiply(xx) as SparseDoubleArray;
            var newF = shift.T.Multiply(F) as SparseDoubleArray;

            var T = newMat.GetSliceDeep(0, L1 * 3 - 1, 0, L1 * 3 - 1);
            var D = newMat.GetSliceDeep(0, L1 * 3 - 1, L1 * 3, _listNode.Count * 3 - 1);
            var fx = newxx.GetSliceDeep(L1 * 3, _listNode.Count * 3 - 1, 0, 0);
            newF = newF.GetSliceDeep(0, L1 * 3 - 1, 0, 0);
            var solve = new SparseLU(T);
            /*var solve = new SparseSVD(T);
            var _U = solve.U;
            var _V = solve.V;
            var _D = solve.D;
            double tol=0.00000001;
            for (int i = 0; i < L1*3; i++)
            {
                if (_D[i, i] > tol) _D[i, i] = 1d / _D[i, i];
                else if (_D[i, i] < -tol) _D[i, i] = 1d / _D[i, i];
                else _D[i, i] = 0d;
            }
            var inv = _V * _D * _U.T;
            */
            var df = D * fx as SparseDoubleArray;
            var b = DoubleArray.From((-newF - df));

            var sol=solve.Solve(b);
            //var sol = inv * b;
            //var sol = _V * _D * _U.T*b;
            var exSol = new SparseDoubleArray(sol.GetLength(0)+fx.GetLength(0),1);
            for (int i = 0; i < L1; i++)
            {
                exSol[i * 3 + 0, 0] = sol[i * 3 + 0, 0];
                exSol[i * 3 + 1, 0] = sol[i * 3 + 1, 0];
                exSol[i * 3 + 2, 0] = sol[i * 3 + 2, 0];
            }
            for (int i = L1; i < _listNode.Count; i++)
            {
                exSol[i * 3 + 0, 0] = fx[(i - L1) * 3 + 0, 0];
                exSol[i * 3 + 1, 0] = fx[(i - L1) * 3 + 1, 0];
                exSol[i * 3 + 2, 0] = fx[(i - L1) * 3 + 2, 0];
            }
            exSol = shift.Multiply(exSol) as SparseDoubleArray;
            foreach (var branch in _listBranch)
            {
                branch.shellCrv = branch.crv.Duplicate() as NurbsCurve;
                for (int i = 0; i < branch.N; i++)
                {
                    var P = branch.shellCrv.Points[i].Location;
                    if (fix)
                    {
                        branch.shellCrv.Points.SetPoint(i, new Point3d(P.X, P.Y, exSol[branch.globalIndex[i] * 3 + 2, 0]));
                    }
                    else
                    {
                        branch.shellCrv.Points.SetPoint(i, new Point3d(exSol[branch.globalIndex[i] * 3 + 0, 0], exSol[branch.globalIndex[i] * 3 + 1, 0], exSol[branch.globalIndex[i] * 3 + 2, 0]));
                    }
                    //if you don't want to allow movements of x and y coordinates, use the following instead of the above. 
                    //branch.shellCrv.Points.SetPoint(i, new Point3d(P.X, P.Y, exSol[branch.globalIndex[i] * 3 + 2, 0]));
                }
            }
            foreach (var leaf in _listLeaf)
            {
                leaf.shellSrf = leaf.srf.Duplicate() as NurbsSurface;
                for (int i = 0; i < leaf.nU; i++)
                {
                    for (int j = 0; j < leaf.nV; j++)
                    {
                        var P = leaf.shellSrf.Points.GetControlPoint(i, j).Location;
                        if (fix)
                        {
                            leaf.shellSrf.Points.SetControlPoint(i, j, new ControlPoint(P.X, P.Y, exSol[leaf.globalIndex[i + j * leaf.nU] * 3 + 2, 0]));
                        }
                        else
                        {
                            leaf.shellSrf.Points.SetControlPoint(i, j, new ControlPoint(exSol[leaf.globalIndex[i + j * leaf.nU] * 3 + 0, 0], exSol[leaf.globalIndex[i + j * leaf.nU] * 3 + 1, 0], exSol[leaf.globalIndex[i + j * leaf.nU] * 3 + 2, 0]));
                        }
                        //if you don't want to allow movements of x and y coordinates, use the following instead of the above. 
                        //leaf.shellSrf.Points.SetControlPoint(i, j, new ControlPoint(P.X, P.Y, exSol[leaf.globalIndex[i + j * leaf.nU] * 3 + 2, 0]));
                    }
                }
            }

        }
        public void mosek1(List<leaf> _listLeaf, List<branch> _listBranch, Dictionary<string, slice> _listSlice,Dictionary<string,range>_listRange,Dictionary<string,range>_listRangeOpen,Dictionary<string,range> _listRangeLeaf, bool obj, double allow)
        {
            // Since the value infinity is never used, we define
            // 'infinity' symbolic purposes only
            double infinity = 0;
            int[] csub = new int[3];// for cones
            int numvar = 0;
            int numcon = 0;
            foreach (var leaf in _listLeaf)
            {
                leaf.varOffset = numvar;
                leaf.conOffset = numcon;
                numvar += (leaf.nU * leaf.nV) + leaf.r * 4;  //z,H11,H22,H12 mean curvature
                numcon += leaf.r * 4;// H11,H22,H12
                if (obj) numvar += leaf.r * 3; //z,target_z, z-_z
                if (obj) numcon += leaf.r * 2; //z, z-target_z
            }

            foreach (var branch in _listBranch)
            {
                branch.varOffset = numvar;
                branch.conOffset = numcon;
                numvar += branch.N + branch.tuples.Count(); //z,D
                if (branch.branchType == branch.type.kink)
                {
                    numcon += 2 * branch.N;//branch->left and right sides
                }
                else if (branch.branchType == branch.type.reinforce||branch.branchType==branch.type.open)
                {
                    numcon += 1 * branch.N; //z=-ax-by-d
                    numcon += 1 * branch.N; //branch->edge(target) 
                }
                else//free
                {
                    numcon += 1 * branch.N; //branch->edge(target)
                }
                numcon += branch.tuples.Count();// D(kink angle)
            }

            foreach (var slice in _listSlice.Values)
            {
                slice.varOffset = numvar;
                slice.conOffset = numcon;
                numvar += 3;  //a,b,d
                if (slice.sliceType == slice.type.fx)
                {
                    numcon++;
                }
            }

            if (obj)
            {
                numvar++;
            }
            //variable settings
            mosek.boundkey[] bkx = new mosek.boundkey[numvar];
            double[] blx = new double[numvar];
            double[] bux = new double[numvar];
            foreach (var leaf in _listLeaf)
            {
                //z
                for (int i = 0; i < leaf.nU * leaf.nV; i++)
                {
                    bkx[i + leaf.varOffset] = mosek.boundkey.fr;
                    blx[i + leaf.varOffset] = -infinity;
                    bux[i + leaf.varOffset] = infinity;
                }
                //H11,H22,H12
                for (int i = 0; i < leaf.r; i++)
                {
                    int n = i * 3 + (leaf.nU * leaf.nV);
                    bkx[n + leaf.varOffset] = mosek.boundkey.fr;
                    blx[n + leaf.varOffset] = -infinity;
                    bux[n + leaf.varOffset] = infinity;
                    bkx[n + 1 + leaf.varOffset] = mosek.boundkey.fr;
                    blx[n + 1 + leaf.varOffset] = -infinity;
                    bux[n + 1 + leaf.varOffset] = infinity;
                    bkx[n + 2 + leaf.varOffset] = mosek.boundkey.fr;
                    blx[n + 2 + leaf.varOffset] = -infinity;
                    bux[n + 2 + leaf.varOffset] = infinity;
                }
                //later mean curvature will be added here
                //
                for (int i = 0; i < leaf.r; i++)
                {
                    int n = i + leaf.r*3+(leaf.nU * leaf.nV);
                    if (leaf.range.rangeType == range.type.lo)
                    {
                        bkx[n + leaf.varOffset] = mosek.boundkey.lo;
                        blx[n + leaf.varOffset] = leaf.range.lb;
                        bux[n + leaf.varOffset] = 0;
                    }
                    else if (leaf.range.rangeType == range.type.up)
                    {
                        bkx[n + leaf.varOffset] = mosek.boundkey.up;
                        blx[n + leaf.varOffset] = 0;
                        bux[n + leaf.varOffset] = leaf.range.ub;
                    }
                    else
                    {
                        bkx[n + leaf.varOffset] = mosek.boundkey.ra;
                        blx[n + leaf.varOffset] = leaf.range.lb;
                        bux[n + leaf.varOffset] = leaf.range.ub;
                    }

                }
                
                ////////////////
                //target z
                if (obj)
                {
                    //z
                    for (int i = 0; i < leaf.r; i++)
                    {
                        bkx[i + (leaf.nU * leaf.nV) + 4 * leaf.r + leaf.varOffset] = mosek.boundkey.fr;
                        blx[i + (leaf.nU * leaf.nV) + 4 * leaf.r + leaf.varOffset] = 0;
                        bux[i + (leaf.nU * leaf.nV) + 4 * leaf.r + leaf.varOffset] = 0;
                    }
                    //target_z
                    for (int i = 0; i < leaf.r; i++)
                    {
                        bkx[i + (leaf.nU * leaf.nV) + 5 * leaf.r + leaf.varOffset] = mosek.boundkey.fx;
                        //reference multiquadric surface
                        blx[i + (leaf.nU * leaf.nV) + 5 * leaf.r + leaf.varOffset] = globalFunc(leaf.tuples[i].x, leaf.tuples[i].y);
                        bux[i + (leaf.nU * leaf.nV) + 5 * leaf.r + leaf.varOffset] = globalFunc(leaf.tuples[i].x, leaf.tuples[i].y);
                    }
                    //z-target_z
                    for (int i = 0; i < leaf.r; i++)
                    {
                        bkx[i + (leaf.nU * leaf.nV) + 6 * leaf.r + leaf.varOffset] = mosek.boundkey.fr;
                        blx[i + (leaf.nU * leaf.nV) + 6 * leaf.r + leaf.varOffset] = 0;
                        bux[i + (leaf.nU * leaf.nV) + 6 * leaf.r + leaf.varOffset] = 0;
                    }
                }
            }
            foreach(var branch in _listBranch)
            {
                if (branch.branchType == branch.type.reinforce )
                {
                    for (int i = 0; i < branch.N; i++)
                    {
                        bkx[i + branch.varOffset] = mosek.boundkey.fr;
                        blx[i + branch.varOffset] = 0;
                        bux[i + branch.varOffset] = 0;
                    }
                    //kink angle parameter
                    for (int i = 0; i < branch.tuples.Count(); i++)
                    {
                        bkx[branch.N + i + branch.varOffset] = mosek.boundkey.lo;
                        blx[branch.N + i + branch.varOffset] = 0.0;
                        bux[branch.N + i + branch.varOffset] = 0;
                    }
                }
                else if (branch.branchType == branch.type.open)
                {
                    for (int i = 0; i < branch.N; i++)
                    {
                        bkx[i + branch.varOffset] = mosek.boundkey.fr;
                        blx[i + branch.varOffset] = 0;
                        bux[i + branch.varOffset] = 0;
                    }
                    //kink angle parameter
                    for (int i = 0; i < branch.tuples.Count(); i++)
                    {
                        if (branch.range.rangeType == range.type.lo)
                        {
                            bkx[branch.N + i + branch.varOffset] = mosek.boundkey.lo;
                            blx[branch.N + i + branch.varOffset] = branch.range.lb;
                            bux[branch.N + i + branch.varOffset] = 0;
                        }
                        else if (branch.range.rangeType == range.type.up)
                        {
                            bkx[branch.N + i + branch.varOffset] = mosek.boundkey.up;
                            blx[branch.N + i + branch.varOffset] = 0;
                            bux[branch.N + i + branch.varOffset] = branch.range.ub;
                        }
                        else
                        {
                            bkx[branch.N + i + branch.varOffset] = mosek.boundkey.ra;
                            blx[branch.N + i + branch.varOffset] = branch.range.lb;
                            bux[branch.N + i + branch.varOffset] = branch.range.ub;
                        }
                        //bkx[branch.N + i + branch.varOffset] = mosek.boundkey.ra;
                        //blx[branch.N + i + branch.varOffset] = 0;
                        //bux[branch.N + i + branch.varOffset] = 0;
                    }
                }
                else if (branch.branchType == branch.type.kink)
                {
                    for (int i = 0; i < branch.N; i++)
                    {
                        bkx[i + branch.varOffset] = mosek.boundkey.fr;
                        blx[i + branch.varOffset] = -infinity;
                        bux[i + branch.varOffset] = infinity;
                    }
                    //kink angle parameter
                    for (int i = 0; i < branch.tuples.Count(); i++)
                    {
                        if (branch.range.rangeType == range.type.lo)
                        {
                            bkx[branch.N + i + branch.varOffset] = mosek.boundkey.lo;
                            blx[branch.N + i + branch.varOffset] = branch.range.lb;
                            bux[branch.N + i + branch.varOffset] = 0;
                        }
                        else if(branch.range.rangeType == range.type.up)
                        {
                            bkx[branch.N + i + branch.varOffset] = mosek.boundkey.up;
                            blx[branch.N + i + branch.varOffset] = 0;
                            bux[branch.N + i + branch.varOffset] = branch.range.ub;
                        }
                        else
                        {
                            bkx[branch.N + i + branch.varOffset] = mosek.boundkey.ra;
                            blx[branch.N + i + branch.varOffset] = branch.range.lb;
                            bux[branch.N + i + branch.varOffset] = branch.range.ub;
                        }
                    }
                }
                else//free
                {
                    for (int i = 0; i < branch.N; i++)
                    {
                        bkx[i + branch.varOffset] = mosek.boundkey.fr;
                        blx[i + branch.varOffset] = -infinity;
                        bux[i + branch.varOffset] = infinity;
                    }
                    //kink angle parameter
                    for (int i = 0; i < branch.tuples.Count(); i++)
                    {
                        bkx[branch.N + i + branch.varOffset] = mosek.boundkey.fr;
                        blx[branch.N + i + branch.varOffset] = -infinity;
                        bux[branch.N + i + branch.varOffset] = infinity;
                    }
                }
            }
            foreach (var slice in _listSlice.Values)
            {
                if (slice.sliceType == slice.type.fx)
                {
                    //add something!
                    bkx[slice.varOffset] = mosek.boundkey.fx;
                    blx[slice.varOffset] = slice.a;
                    bux[slice.varOffset] = slice.a;
                    bkx[slice.varOffset + 1] = mosek.boundkey.fx;
                    blx[slice.varOffset + 1] = slice.b;
                    bux[slice.varOffset + 1] = slice.b;
                    bkx[slice.varOffset + 2] = mosek.boundkey.fx;
                    blx[slice.varOffset + 2] = slice.d;
                    bux[slice.varOffset + 2] = slice.d;
                }
                else
                {
                    bkx[slice.varOffset] = mosek.boundkey.fr;
                    blx[slice.varOffset] = -infinity;
                    bux[slice.varOffset] = infinity;
                    bkx[slice.varOffset + 1] = mosek.boundkey.fr;
                    blx[slice.varOffset + 1] = -infinity;
                    bux[slice.varOffset + 1] = infinity;
                    bkx[slice.varOffset + 2] = mosek.boundkey.fr;
                    blx[slice.varOffset + 2] = -infinity;
                    bux[slice.varOffset + 2] = infinity;
                }
            }
            if (obj)
            {
                bkx[numvar - 1] = mosek.boundkey.fx;
                blx[numvar - 1] = allow;
                bux[numvar - 1] = allow;
                
                //bkx[numvar - 1] = mosek.boundkey.fr;
                //blx[numvar - 1] = -infinity;
                //bux[numvar - 1] = infinity;
            }

            // Make mosek environment.
            using (mosek.Env env = new mosek.Env())
            {
                // Create a task object.
                using (mosek.Task task = new mosek.Task(env, 0, 0))
                {
                    // Directs the log task stream to the user specified
                    // method msgclass.streamCB
                    task.set_Stream(mosek.streamtype.log, new msgclass(""));

                    /* Append 'numcon' empty constraints.
                       The constraints will initially have no bounds. */
                    task.appendcons(numcon);

                    /* Append 'numvar' variables.
                       The variables will initially be fixed at zero (x=0). */
                    task.appendvars(numvar);

                    for (int j = 0; j < numvar; ++j)
                    {
                        task.putvarbound(j, bkx[j], blx[j], bux[j]);
                    }
                    double root2 = Math.Sqrt(2);
                    foreach (var leaf in listLeaf)
                    {

                        double[] grad = new double[leaf.tuples[0].nNode];
                        double[] grad0 = new double[leaf.tuples[0].nNode];
                        double[] grad1i = new double[leaf.tuples[0].nNode];
                        double[] grad1j = new double[leaf.tuples[0].nNode];
                        //define H11,H12,H22
                        for (int i = 0; i < leaf.r; i++)
                        {
                            int N11 = i * 3; //condition number
                            int N22 = i * 3 + 1;
                            int N12 = i * 3 + 2;
                            int target = i * 3 + (leaf.nU * leaf.nV) + leaf.varOffset;   //variable number
                            task.putaij(N11+leaf.conOffset, target, -1);
                            task.putconbound(N11 + leaf.conOffset, mosek.boundkey.fx, 0, 0);
                            task.putaij(N22 + leaf.conOffset, target + 1, -1);
                            task.putconbound(N22 + leaf.conOffset, mosek.boundkey.fx, 0, 0);
                            task.putaij(N12 + leaf.conOffset, target + 2, -1);
                            task.putconbound(N12 + leaf.conOffset, mosek.boundkey.fx, 0, 0);
                            //N11
                            leaf.tuples[i].d2[0, 0].CopyTo(grad, 0);
                            leaf.tuples[i].d0.CopyTo(grad0, 0);
                            leaf.tuples[i].d1[0].CopyTo(grad1i, 0);
                            leaf.tuples[i].d1[0].CopyTo(grad1j, 0);
                            for (int k = 0; k < leaf.tuples[i].nNode; k++)
                            {
                                for (int j = 0; j < leaf.tuples[i].elemDim; j++)
                                {
                                    grad[k] -= leaf.tuples[i].Gammaijk[0, 0, j] * leaf.tuples[i].d1[j][k];
                                }
                                double val = 0;
                                val += grad[k];
                                task.putaij(N11 + leaf.conOffset, leaf.tuples[i].internalIndex[k] + leaf.varOffset, -val / root2);
                            }
                            //N22
                            leaf.tuples[i].d2[1, 1].CopyTo(grad, 0);
                            leaf.tuples[i].d0.CopyTo(grad0, 0);
                            leaf.tuples[i].d1[1].CopyTo(grad1i, 0);
                            leaf.tuples[i].d1[1].CopyTo(grad1j, 0);
                            for (int k = 0; k < leaf.tuples[i].nNode; k++)
                            {
                                for (int j = 0; j < leaf.tuples[i].elemDim; j++)
                                {
                                    grad[k] -= leaf.tuples[i].Gammaijk[1, 1, j] * leaf.tuples[i].d1[j][k];
                                }
                                double val = 0;
                                val += grad[k];
                                task.putaij(N22 + leaf.conOffset, leaf.tuples[i].internalIndex[k] + leaf.varOffset, -val / root2);
                            }
                            //N12
                            leaf.tuples[i].d2[0, 1].CopyTo(grad, 0);
                            leaf.tuples[i].d0.CopyTo(grad0, 0);
                            leaf.tuples[i].d1[0].CopyTo(grad1i, 0);
                            leaf.tuples[i].d1[1].CopyTo(grad1j, 0);
                            for (int k = 0; k < leaf.tuples[i].nNode; k++)
                            {
                                for (int j = 0; j < leaf.tuples[i].elemDim; j++)
                                {
                                    grad[k] -= leaf.tuples[i].Gammaijk[0, 1, j] * leaf.tuples[i].d1[j][k];
                                }
                                double val = 0;
                                val += grad[k];
                                task.putaij(N12 + leaf.conOffset, leaf.tuples[i].internalIndex[k] + leaf.varOffset, -val);
                            }
                        }
                        // mean curvature will be added here
                        //
                        //
                        for (int i = 0; i < leaf.r; i++)
                        {
                            int target = i * 3 + (leaf.nU * leaf.nV) + leaf.varOffset;   //variable number
                            int target2 = i + leaf.r*3+(leaf.nU * leaf.nV) + leaf.varOffset;   //variable number

                            int NH= leaf.r*3+i; //condition number
                            task.putaij(NH + leaf.conOffset, target, 1);
                            task.putaij(NH + leaf.conOffset, target + 1, 1);
                            task.putaij(NH + leaf.conOffset, target2, -1);
                            task.putconbound(NH+leaf.conOffset, mosek.boundkey.fx, 0, 0);
                         
                        }

                        //if (leaf.leafType == leaf.type.convex)
                        for (int i = 0; i < leaf.r; i++)
                        {
                            int N11 = i * 3 + (leaf.nU * leaf.nV); //variable number
                            int N22 = i * 3 + 1 + (leaf.nU * leaf.nV);
                            int N12 = i * 3 + 2 + (leaf.nU * leaf.nV);

                            csub[0] = N11 + leaf.varOffset;
                            csub[1] = N22 + leaf.varOffset;
                            csub[2] = N12 + leaf.varOffset;
                            task.appendcone(mosek.conetype.rquad,
                                            0.0, // For future use only, can be set to 0.0 
                                            csub);
                            /*if (obj2)
                            {
                                task.putcj(N11, leaf.tuples[i].Gij[0, 0]);
                                task.putcj(N22, leaf.tuples[i].Gij[1, 1]);
                                task.putcj(N12, 2*leaf.tuples[i].Gij[0, 1]);
                            }*/
                        }
                        if (obj)
                        {
                            double[] grad00 = new double[leaf.tuples[0].nNode];
                            for (int i = 0; i < leaf.r; i++)
                            {
                                leaf.tuples[i].d0.CopyTo(grad00, 0);
                                for (int k = 0; k < leaf.tuples[i].nNode; k++)
                                {
                                    task.putaij(leaf.conOffset + leaf.r * 4 + i, leaf.varOffset + leaf.tuples[i].internalIndex[k], grad00[k]);
                                }
                                task.putaij(leaf.conOffset + leaf.r * 4 + i, leaf.varOffset + leaf.nU*leaf.nV+leaf.r*4+i,-1);
                                task.putconbound(leaf.conOffset + leaf.r * 4 + i, mosek.boundkey.fx, 0, 0);
                            }
                            for (int i = 0; i < leaf.tuples.Count(); i++)
                            {
                                task.putaij(leaf.conOffset + leaf.r * 5 + i, leaf.varOffset + leaf.nU * leaf.nV + leaf.r * 4 + i, 1);
                                task.putaij(leaf.conOffset + leaf.r * 5 + i, leaf.varOffset + leaf.nU * leaf.nV + leaf.r * 5 + i, -1);
                                task.putaij(leaf.conOffset + leaf.r * 5 + i, leaf.varOffset + leaf.nU * leaf.nV + leaf.r * 6 + i, -1);
                                task.putconbound(leaf.conOffset + leaf.r * 5 + i, mosek.boundkey.fx, 0, 0);
                            }
                        }
                    }
                    
                    if (obj)
                    {
                        List<int> dsub=new List<int>();
                        dsub.Add(numvar-1);
                        foreach (var leaf in _listLeaf)
                        {
                            for (int i = 0; i < leaf.r; i++)
                            {
                                dsub.Add(leaf.varOffset + leaf.nU * leaf.nV + leaf.r * 6 + i);
                            }
                        }
                        task.appendcone(mosek.conetype.quad, 0.0, dsub.ToArray());
                    }
                    foreach (var branch in _listBranch)
                    {
                        if (branch.branchType == branch.type.kink)
                        {
                            tieBranchD1(branch, branch.left, task, 2, 0);
                            tieBranchD1(branch, branch.right, task, 2, 1);
                            defineKinkAngle2(branch,branch.left,branch.right,task, branch.conOffset + branch.N*2, branch.varOffset + branch.N);
                            /*if (branch.obj)
                            {
                                for (int i = 0; i < branch.tuples.Count(); i++)
                                {
                                    task.putcj(branch.N + i + branch.varOffset, 1);
                                }
                            }*/
                        }
                        else if (branch.branchType == branch.type.reinforce || branch.branchType == branch.type.open)
                        {
                            int iA = _listSlice[branch.sliceKey].varOffset;
                            int iB = _listSlice[branch.sliceKey].varOffset + 1;
                            int iD = _listSlice[branch.sliceKey].varOffset + 2;
                            //height parameter
                            for (int i = 0; i < branch.N; i++)
                            {
                                double x = branch.crv.Points[i].Location.X;
                                double y = branch.crv.Points[i].Location.Y;
                                task.putconbound(branch.conOffset + branch.N + branch.tuples.Count() + i, mosek.boundkey.fx, 0, 0);
                                task.putaij(branch.conOffset + branch.N + branch.tuples.Count() + i, branch.varOffset + i, 1);//z
                                task.putaij(branch.conOffset + branch.N + branch.tuples.Count() + i, iA, x);//ax
                                task.putaij(branch.conOffset + branch.N + branch.tuples.Count() + i, iB, y);//by
                                task.putaij(branch.conOffset + branch.N + branch.tuples.Count() + i, iD, 1);//d
                            }
                            tieBranchD1(branch, branch.target, task, 1, 0);
                            defineKinkAngleC(branch, branch.target, task, branch.conOffset + branch.N, branch.varOffset + branch.N);
                        }
                        else
                        {
                            tieBranchD1(branch, branch.target, task, 1, 0);
                            defineKinkAngle(branch,branch.target, task, branch.conOffset + branch.N, branch.varOffset + branch.N);
                        }
                    }
                    //task.putcj(numvar - 1, 1);
                    task.putintparam(mosek.iparam.intpnt_max_iterations, 200000000);//20000000
                    task.putintparam(mosek.iparam.intpnt_solve_form, mosek.solveform.dual);
                    task.putobjsense(mosek.objsense.minimize);
                    //task.writedata("c:/out/mosek_task_dump.opf");


                    task.optimize();
                    // Print a summary containing information
                    //   about the solution for debugging purposes
                    task.solutionsummary(mosek.streamtype.msg);

                    mosek.solsta solsta;
                    /* Get status information about the solution */
                    task.getsolsta(mosek.soltype.itr, out solsta);

                    double[] xx = new double[numvar];

                    task.getxx(mosek.soltype.itr, // Basic solution.     
                                    xx);

                    switch (solsta)
                    {
                        case mosek.solsta.optimal:
                            System.Windows.Forms.MessageBox.Show("Optimal primal solution\n");
                            break;
                        case mosek.solsta.near_optimal:
                            System.Windows.Forms.MessageBox.Show("Near Optimal primal solution\n");
                            break;
                        case mosek.solsta.dual_infeas_cer:
                        case mosek.solsta.prim_infeas_cer:
                        case mosek.solsta.near_dual_infeas_cer:
                        case mosek.solsta.near_prim_infeas_cer:
                            Console.WriteLine("Primal or dual infeasibility.\n");
                            break;
                        case mosek.solsta.unknown:
                            System.Windows.Forms.MessageBox.Show("Unknown solution status\n");
                            break;
                        default:
                            System.Windows.Forms.MessageBox.Show("Other solution status\n");
                            break;

                    }
                    //store airy potential
                    System.Windows.Forms.MessageBox.Show(string.Format("error={0}", xx[numvar - 1]));
                    foreach (var leaf in listLeaf)
                    {
                        double[] x = new double[leaf.nU * leaf.nV];
                        for (int j = 0; j < leaf.nV; j++)
                        {
                            for (int i = 0; i < leaf.nU; i++)
                            {
                                x[i + j * leaf.nU] = xx[i + j * leaf.nU + leaf.varOffset];
                            }
                        }
                        leaf.myMasonry.setupAiryPotentialFromList(x);
                    }
                    foreach (var leaf in listLeaf)
                    {
                        foreach (var tup in leaf.tuples)
                        {
                            leaf.myMasonry.elemList[tup.index].computeStressFunction(tup);
                        }
                    }
                    foreach (var branch in _listBranch)
                    {
                        double[] x = new double[branch.N];
                        for (int i = 0; i < branch.N; i++)
                        {
                            x[i] = xx[i + branch.varOffset];
                        }
                        branch.myArch.setupAiryPotentialFromList(x);
                    }
                    foreach (var slice in _listSlice.Values)
                    {
                        slice.a = xx[slice.varOffset];
                        slice.b = xx[slice.varOffset + 1];
                        slice.d = xx[slice.varOffset + 2];
                        double norm = Math.Sqrt(slice.a * slice.a + slice.b * slice.b + 1);
                        var pl = new Rhino.Geometry.Plane(slice.a, slice.b, 1d, slice.d / norm);
                        slice.update(pl);
                    }
                    foreach (var branch in listBranch)
                    {
                        foreach (var tup in branch.tuples)
                        {
                            if (branch.branchType == branch.type.kink)
                            {
                                branch.left.myMasonry.elemList[tup.left.index].computeTangent(tup.left);
                                branch.right.myMasonry.elemList[tup.right.index].computeTangent(tup.right);
                            }
                            else if (branch.branchType == branch.type.fix)
                            {
                                branch.target.myMasonry.elemList[tup.target.index].computeTangent(tup.target);
                            }
                            else
                            {

                                branch.target.myMasonry.elemList[tup.target.index].computeTangent(tup.target);
                                var vars = branch.slice.pl.GetPlaneEquation();
                                branch.target.myMasonry.elemList[tup.target.index].computeTangent(tup.target, vars[0], vars[1], vars[2], vars[3]); //valDc
                            }
                        }
                    }
                    foreach (var leaf in _listLeaf)
                    {
                        for (int i = 0; i < leaf.r; i++)
                        {
                            int target = i + leaf.r*3+(leaf.nU * leaf.nV) + leaf.varOffset;
                            leaf.tuples[i].NH = xx[target];
                        }
                    }
                    foreach (var branch in _listBranch)
                    {
                        branch.airyCrv = branch.crv.Duplicate() as NurbsCurve;
                        for (int j = 0; j < branch.N; j++)
                        {
                            var P = branch.crv.Points[j];
                            branch.airyCrv.Points.SetPoint(j, new Point3d(P.Location.X, P.Location.Y, xx[j + branch.varOffset]));
                        }
                        for (int i = 0; i < branch.tuples.Count(); i++)
                        {
                            //branch.tuples[i].z = branch.airyCrv.PointAt(branch.tuples[i].t).Z;
                            //int D = i + branch.N;
                            if (branch.branchType == branch.type.open)
                            {
                                branch.tuples[i].H[0, 0] = branch.tuples[i].target.valD - branch.tuples[i].target.valDc;
                            }
                            else if (branch.branchType == branch.type.reinforce)
                            {
                                branch.tuples[i].H[0, 0] = branch.tuples[i].target.valD - branch.tuples[i].target.valDc;
                            }
                            else if (branch.branchType == branch.type.fix)
                            {
                                branch.tuples[i].H[0, 0] = 0;
                            }
                            else
                            {
                                //afterwards, check why these two values do not match.
                                //branch.tuples[i].H[0, 0] = branch.tuples[i].left.valD + branch.tuples[i].right.valD;
                                branch.tuples[i].H[0, 0] = xx[branch.N + i + branch.varOffset];
                            }
                        }
                    }
                    foreach (var range in _listRangeLeaf.Values)
                    {
                        double min = 10000d, max = -10000d;
                        foreach (var leaf in range.lL)
                        {
                            for (int i = 0; i < leaf.tuples.Count(); i++)
                            {
                                if (leaf.tuples[i].NH > max) max = leaf.tuples[i].NH;
                                if (leaf.tuples[i].NH < min) min = leaf.tuples[i].NH;
                            }
                        }
                        range.lastMin = min;
                        range.lastMax = max;
                        range.firstPathDone = true;
                    }
                    foreach (var range in _listRange.Values)
                    {
                        double min=10000d,max=-10000d;
                        foreach (var branch in range.lB)
                        {
                            for(int i=0;i<branch.tuples.Count();i++)
                            {
                                if (branch.tuples[i].H[0, 0] > max) max = branch.tuples[i].H[0, 0];
                                if (branch.tuples[i].H[0, 0] < min) min = branch.tuples[i].H[0, 0];
                            }
                        }
                        range.lastMin = min;
                        range.lastMax = max;
                        range.firstPathDone = true;
                    }
                    foreach (var range in _listRangeOpen.Values)
                    {
                        double min = 10000d, max = -10000d;
                        foreach (var branch in range.lB)
                        {
                            for (int i = 0; i < branch.tuples.Count(); i++)
                            {
                                if (branch.tuples[i].H[0, 0] > max) max = branch.tuples[i].H[0, 0];
                                if (branch.tuples[i].H[0, 0] < min) min = branch.tuples[i].H[0, 0];
                            }
                        }
                        range.lastMin = min;
                        range.lastMax = max;
                        range.firstPathDone = true;
                    }

                    foreach (var leaf in _listLeaf)
                    {
                        leaf.airySrf = leaf.srf.Duplicate() as NurbsSurface;
                        for (int j = 0; j < leaf.nV; j++)
                        {
                            for (int i = 0; i < leaf.nU; i++)
                            {
                                var P = leaf.srf.Points.GetControlPoint(i, j);
                                leaf.airySrf.Points.SetControlPoint(i, j, new ControlPoint(P.Location.X, P.Location.Y, xx[i + j * leaf.nU + leaf.varOffset]));
                            }
                        }
                    }
                }
            }
        }
        void hodgeStar(List<leaf> _listLeaf, List<branch> _listBranch, Func<double, double> coeff,double sScale)
        {
            foreach (var branch in _listBranch)
            {
                for (int i = 0; i < branch.tuples.Count(); i++)
                {
                    double g = branch.tuples[i].gij[0, 0];
                    double val = coeff(g);
                    branch.tuples[i].SPK[0, 0] = branch.tuples[i].H[0, 0] * val*sScale;
                    if(branch.tuples[i].SPK[0,0]<0)
                    {
                        branch.tuples[i].SPK[0, 0] = 0d;
                    }
                }
            }
            foreach (var leaf in _listLeaf)
            {
                for (int j = 0; j < leaf.r; j++)
                {
                    //Hodge star
                    double g = leaf.tuples[j].refDv * leaf.tuples[j].refDv;

                    leaf.tuples[j].SPK[0, 0] = leaf.tuples[j].H[1, 1] / g;
                    leaf.tuples[j].SPK[1, 1] = leaf.tuples[j].H[0, 0] / g;
                    leaf.tuples[j].SPK[0, 1] = -leaf.tuples[j].H[0, 1] / g;
                    leaf.tuples[j].SPK[1, 0] = -leaf.tuples[j].H[1, 0] / g;
                    leaf.tuples[j].computeEigenVectors();
                    var tup = leaf.tuples[j];
                    var det = tup.SPK[0, 0] * tup.SPK[1, 1] - tup.SPK[0, 1] * tup.SPK[1, 0];
                    if (tup.eigenValues[0] < 0 || tup.eigenValues[1] < 0)
                    {
                        if (tup.eigenValues[0] < 0) tup.eigenValues[0] = 0;
                        if (tup.eigenValues[1] < 0) tup.eigenValues[1] = 0;
                        //P
                        double A11 = tup.eigenVectorsB[0][0];
                        double A12 = tup.eigenVectorsB[0][1];
                        double A21 = tup.eigenVectorsB[1][0];
                        double A22 = tup.eigenVectorsB[1][1];
                        double det2 = A11 * A22 - A12 * A21;
                        //P^-1
                        double B11 = A22 / det2;
                        double B22 = A11 / det2;
                        double B12 = -A12 / det2;
                        double B21 = -A21 / det2;
                        double C11 = B11 * tup.eigenValues[0];
                        double C12 = B12 * tup.eigenValues[1];
                        double C21 = B21 * tup.eigenValues[0];
                        double C22 = B22 * tup.eigenValues[1];
                        double D11 = C11 * A11 + C12 * A21;
                        double D12 = C11 * A12 + C12 * A22;
                        double D21 = C21 * A11 + C22 * A21;
                        double D22 = C21 * A12 + C22 * A22;
                        tup.SPK[0, 0] = D11 * tup.Gij[0, 0] + D12 * tup.Gij[1, 0];
                        tup.SPK[0, 1] = D11 * tup.Gij[0, 1] + D12 * tup.Gij[1, 1];
                        tup.SPK[1, 0] = D21 * tup.Gij[0, 0] + D22 * tup.Gij[1, 0];
                        tup.SPK[1, 1] = D21 * tup.Gij[1, 0] + D22 * tup.Gij[1, 1];
                    }
                    tup.SPK[0, 0] *= sScale;
                    tup.SPK[1, 0] *= sScale;
                    tup.SPK[0, 1] *= sScale;
                    tup.SPK[1, 1] *= sScale;
                }
            }
            //For visualization
            crossMagenta.Clear();
            crossCyan.Clear();
            foreach (var leaf in listLeaf)
            {
                foreach (var tuple in leaf.tuples)
                {
                    for (int i = 0; i < 2; i++)
                    {
                        if (tuple.eigenValues[i] < 0)
                        {
                            double s = tuple.eigenValues[i]*sScale;
                            //double s = 0.1;
                            Point3d S = new Point3d(tuple.x - tuple.eigenVectors[i][0] * s, tuple.y - tuple.eigenVectors[i][1] * s, tuple.z - tuple.eigenVectors[i][2] * s);
                            Point3d E = new Point3d(tuple.x + tuple.eigenVectors[i][0] * s, tuple.y + tuple.eigenVectors[i][1] * s, tuple.z + tuple.eigenVectors[i][2] * s);
                            Line line = new Line(S, E);
                            line.Transform(zDown);
                            crossCyan.Add(line);
                        }
                        else
                        {
                            double s = tuple.eigenValues[i]*sScale;
                            //double s = 0.1;
                            Point3d S = new Point3d(tuple.x - tuple.eigenVectors[i][0] * s, tuple.y - tuple.eigenVectors[i][1] * s, tuple.z - tuple.eigenVectors[i][2] * s);
                            Point3d E = new Point3d(tuple.x + tuple.eigenVectors[i][0] * s, tuple.y + tuple.eigenVectors[i][1] * s, tuple.z + tuple.eigenVectors[i][2] * s);
                            Line line = new Line(S, E);
                            line.Transform(zDown);
                            crossMagenta.Add(line);
                        }
                    }
                }
            }

        }
    }
}
