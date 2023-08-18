from IPython.display import HTML


def hide_code():
    """
    This function hides the code cells in a jupyter notebook. It helps to a
    cleaner in visualization of results.

    Note: this requires
    ipython
    ipywidgets
    """
    html = HTML('''<script>
    code_show=true;
    function code_toggle() {
    if (code_show){
    $('div.input').hide();
    } else {
    $('div.input').show();
    }
    code_show = !code_show
    }
    $( document ).ready(code_toggle);
    </script>
    <form action="javascript:code_toggle()"><input type="submit"
    value="Click here to toggle on/off the raw code."></form>''')
    return html


"""
tips:

to have interactive plots:
   %matplotlib notebook

to remove outputs completely:
   %%capture

to run a bash line:
   ! <command>
"""
