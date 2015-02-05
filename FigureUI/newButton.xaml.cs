using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace Mothra.UI
{
    /// <summary>
    /// newButton.xaml の相互作用ロジック
    /// </summary>
    public partial class newButton : UserControl
    {
        public Action function=null;
        public newButton(Action func)
        {
            function = func;
            InitializeComponent();
        }
        public newButton()
        {
            InitializeComponent();
        }
        public string Text
        {
            set
            {
                this.button1.Content = value;
            }
            get
            {
                return (string)this.button1.Content;
            }
        }
        private void button_Click(object sender, RoutedEventArgs e)
        {
            if(function!=null)function();
        }
    }
}
