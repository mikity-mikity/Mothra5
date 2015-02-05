using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using mikity;

namespace Simplex
{
    class Program
    {
        [STAThread]
        static void Main(string[] args)
        {
            int hSeg = 16;
            int vSeg = 8;
            mikity.index myIndex = new mikity.index();                  //節点番号と座標変数の対応表
            mikity.index mask = new mikity.index();                  //節点番号と座標変数の対応表
            mikity.shape myShape = new mikity.shape();                  //座標変数
            mikity.elements Triangles = new mikity.elements();          //三角形要素
            mikity.elements myElements = new mikity.elements();         //張力要素リスト
            mikity.elements Border = new mikity.elements();         //境界要素リスト

            //形状定義
            System.Random rand = new System.Random(0);
            for (int i = 0; i < hSeg; i++)
            {
                for (int j = 0; j < vSeg + 1; j++)
                {
                    int num = i + j * (hSeg);
                    myIndex.Add(new int[3] { num * 3, num * 3 + 1, num * 3 + 2 });
                    myShape.AddRange(new double[3] { (rand.NextDouble() - 0.5) * 15d, (rand.NextDouble() - 0.5) * 15d, (rand.NextDouble() - 0.5) * 15d });
                }
            }
            double Radius = 10.0;
            for (int i = 0; i < hSeg; i++)
            {
                int num = i + 0 * (hSeg);
                mask.Add(myIndex[num]);
                double theta = System.Math.PI * 2 / hSeg * i;
                myShape[mask.Last()] = new double[3] { Radius*System.Math.Cos(theta), Radius*System.Math.Sin(theta), -2.0 };
            }
            for (int i = 0; i < hSeg; i++)
            {
                int num = i + vSeg * (hSeg);
                double theta = System.Math.PI * 2 / hSeg * i;
                mask.Add(myIndex[num]);
                myShape[mask.Last()] = new double[3] { Radius * System.Math.Cos(theta), Radius * System.Math.Sin(theta), 0.0 };
            }
            //基準正方形を構成する2つの三角形
            int[,] _t = (matrix_INT)(new int[2, 2] { { 0, 1 }, { hSeg, hSeg+1 } });
            Triangles.Add(new int[3] { _t[0, 0], _t[0, 1], _t[1, 0] });
            Triangles.Add(new int[3] { _t[0, 1], _t[1, 1], _t[1, 0] });
            int[] tmp = vector.lin(0, 1);
            //作成した基準正方形を並べて行を作成
            for (int i = 1; i < hSeg; i++)
            {
                Triangles.AddRange((elements)Triangles[tmp] + i);
            }
            Triangles[Triangles.Count-2][1] = 0;
            Triangles[Triangles.Count-1][0] = 0;
            Triangles[Triangles.Count-1][1] = hSeg;
            //作成した基準行を並べて膜を作成
            tmp = vector.lin(0, 2*hSeg-1);
            for (int i = 1; i < vSeg; i++)
            {
                Triangles.AddRange((elements)Triangles[tmp] + i * (hSeg));

            }
            myElements.AddRange(Triangles);

            //3項法の準備
            mikity.term p = new mikity.term(myShape.Count);          //加速度
            mikity.term q = new mikity.term(myShape.Count);          //速度
            mikity.term x = new mikity.term(myShape);                //位置

            //GUIと可変パラメータの準備
            mikity.visualize.FigureUI figUI = new mikity.visualize.FigureUI();
            //時間の刻み幅
            mikity.visualize.slider s1 = figUI.addSlider(min: 1, step: 2, max: 44, val: 38, text: "dt");
            s1.Converter = val => Math.Pow(10, ((double)val / 10) - 4.8) * 2.0;
            System.Collections.Generic.List<double> dt = new System.Collections.Generic.List<double>();
            //リング間の距離
            mikity.visualize.slider s2 = figUI.addSlider(min: 0, step: 2, max: 4000, val: 800, text: "Ring");
            s2.Converter = val => (val)/200;
            //描画ウインドウの準備
            mikity.visualize.figure3D myFigure = new mikity.visualize.figure3D();
            myFigure.xtitle = "Tanzbrunnen";
            System.Collections.Generic.List<double> fps = new System.Collections.Generic.List<double>();
            System.Collections.Generic.List<int> ticks = new System.Collections.Generic.List<int>();
            ticks.AddRange(new int[101]);
            //3項法のメインループ
            int t = 0;
            mikity.visualize.drawelements.scale(0.5);
            while (figUI.Exit == false)
            {
                if (figUI.Pause == false)
                {
                    t++;
                    ticks.Add(System.Environment.TickCount);
                    ticks.Remove(ticks[0]);
                    fps.Add((double)1000 / (ticks[100] - ticks[0]) * 100);
                    dt.Add(s1.value);
                    double w1 = s2.value;//リング
                    for (int i = 0; i < hSeg; i++)
                    {
                        double theta = System.Math.PI * 2 / hSeg * i;
                        myShape[mask[i]] = new double[3] { Radius * System.Math.Cos(theta), Radius * System.Math.Sin(theta), -w1 };
                    }
                    for (int i = 0; i < hSeg; i++)
                    {
                        double theta = System.Math.PI * 2 / hSeg * i;
                        myShape[mask[i + hSeg]] = new double[3] { Radius * System.Math.Cos(theta), Radius * System.Math.Sin(theta), -0.0 };
                    }
                    if (figUI.Exit == true) break;
                    System.Threading.Tasks.Parallel.ForEach(myElements, (element e) =>
                    {
                        e.Update(myShape, myIndex);
                        e.Stress = e.Metric;
                        e.makeOmega();
                    });
                    vector grad = myElements.Gradient();
                    vector.mask(grad, mask);
                    p.Add(grad);
                    if (vector.norm(p[t]) != 0)
                    {
                        p[t] = (vector)p[t] / vector.norm(p[t]);
                    }
                    q.Add((0.98 * (vector)q[t - 1]) - (vector)p[t] * dt[t - 1]);
                    //x.Add(myShape + (vector)q[t] * dt[t - 1]);
                    x.Add(myShape - (vector)p[t] * dt[t - 1]*10);
                    myShape = x[t];
                    if (t % 5 == 0)
                    {
                        myFigure.drawlater();
                        myFigure.clf();
                        myFigure.WriteLine("dt=" + dt[t - 1].ToString("0.000e00"));
                        myFigure.WriteLine("FPS=" + fps[t - 1].ToString("0.0"));
                        myFigure.set_shape(myShape, myIndex);
                        mikity.visualize.drawelements.draw_fixedpoints(myFigure, mask, myIndex);
                        mikity.visualize.drawelements.draw_triangles(myFigure, Triangles, myIndex);
                        myFigure.drawnow();
                    }
                }
                if (figUI.expPNG == true)
                {
                    myFigure.exportPNG();
                }
                if (figUI.expEMF == true)
                {
                    myFigure.exportEMF();
                }
                if (figUI.expSeq == true)
                {
                    System.Windows.MessageBoxResult res = System.Windows.MessageBox.Show("Click OK to export sequential BMPs", "Message from mikity", System.Windows.MessageBoxButton.OKCancel, System.Windows.MessageBoxImage.Information);
                    if (res == System.Windows.MessageBoxResult.OK)
                    {
                        mikity.visualize.sequence.init();
                        using (mikity.visualize.figure3D myFigure2 = new mikity.visualize.figure3D())
                        {
                            myFigure2.xtitle = "Tanzbrunnen";
                            for (int i = 1; i <= t; i += 5)
                            {
                                myFigure2.drawlater();
                                myFigure2.clf();
                                myFigure2.WriteLine("dt=" + dt[i - 1].ToString("0.000e00"));
                                myFigure2.WriteLine("FPS=" + fps[i - 1].ToString("0.0"));
                                myFigure2.set_shape(x[i], myIndex);
                                mikity.visualize.drawelements.draw_fixedpoints(myFigure2, mask, myIndex);
                                mikity.visualize.drawelements.draw_triangles(myFigure2, Triangles, myIndex);
                                myFigure2.drawnow();
                                myFigure2.exportBMP((int)i / 5);
                                System.Windows.Forms.Application.DoEvents();
                            }
                        }
                    }
                    System.Windows.MessageBox.Show("Finish!", "Message from mikity", System.Windows.MessageBoxButton.OK, System.Windows.MessageBoxImage.Information);
                }
                System.Windows.Forms.Application.DoEvents();
            }


        }
    }
}
