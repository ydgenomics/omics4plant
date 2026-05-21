# 导入 Dash 相关库
import dash
from dash import dcc, html  # dcc 是 Dash 核心组件库，html 是 HTML 组件库

# 创建一个 Dash 应用实例
app = dash.Dash(__name__)

# 定义应用的布局
app.layout = html.Div(children=[
    # 添加一个标题
    html.H1(children='你好，Dash！'),

    # 添加一段描述文字
    html.Div(children='''
        Dash：一个用于 Python 的 Web 应用框架。
    '''),

    # 添加一个图表
    dcc.Graph(
        id='example-graph',  # 图表的 ID，用于回调函数
        figure={
            'data': [  # 图表的数据
                {'x': [1, 2, 3], 'y': [4, 1, 2], 'type': 'bar', 'name': '上海'},
                {'x': [1, 2, 3], 'y': [2, 4, 5], 'type': 'bar', 'name': '北京'},
            ],
            'layout': {  # 图表的布局
                'title': 'Dash 数据可视化示例'  # 图表的标题
            }
        }
    )
])

# 运行应用
if __name__ == '__main__':
    app.run(debug=True)  # 启动应用，debug=True 表示开启调试模式