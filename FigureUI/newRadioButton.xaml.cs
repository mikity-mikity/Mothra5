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
    public partial class newRadioButton : UserControl
    {
        public Action<int> function=null;
        public int number = -1;
        public newRadioButton(Action<int> func)
        {
            function = func;
            InitializeComponent();
        }
        public newRadioButton()
        {
            InitializeComponent();
        }
        public string Text
        {
            set
            {
                this.radiobutton1.Content = value;
            }
            get
            {
                return (string)this.radiobutton1.Content;
            }
        }
 
        private void button1_Checked(object sender, RoutedEventArgs e)
        {
            if (function != null) function(this.number);
        }
    }
}
