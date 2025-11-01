import gradio as gr
import pandas as pd
from pridict2_pegRNA_design import findpegRNAs

def run_pridict(sequence, editor, pbs_min, pbs_max, rt_min, rt_max):
    try:
        # Create a dummy DataFrame row as input for the function
        df_row = pd.Series({
            'editseq': sequence,
            'sequence_name': 'gradio_input'
        })
        
        # Define the ranges for PBS and RT lengths
        pbs_range = list(range(int(pbs_min), int(pbs_max) + 1))
        rt_range = list(range(int(rt_min), int(rt_max) + 1))

        # Run the main function from your script
        df_pegrna, df_nicking, error_message = findpegRNAs(
            df_row,
            PBSlengthrange=pbs_range,
            RToverhanglengthrange=rt_range,
            editor=editor
        )

        if error_message:
            return None, None, f"An error occurred: {error_message}"

        # Return the results as Gradio DataFrames and a status message
        return df_pegrna, df_nicking, "Successfully generated pegRNAs!"

    except Exception as e:
        return None, None, f"A critical error occurred: {str(e)}"

# --- Gradio Interface Definition ---
with gr.Blocks() as demo:
    gr.Markdown("# PRIDICT2 pegRNA Designer")
    gr.Markdown("Run the modified PRIDICT2 script to design pegRNAs with custom PAM support.")

    with gr.Row():
        editor_dropdown = gr.Dropdown(
            label="Select Prime Editor (and PAM sequence)",
            choices=['PE2-NGG', 'PE2-NG', 'PE2-SpRY', 'PE2-NRN'],
            value='PE2-NRN'
        )
    
    sequence_input = gr.Textbox(
        label="Input Sequence",
        placeholder="Enter your sequence with the edit in brackets, e.g., ...GCA[G/A]TGC..."
    )
    
    with gr.Accordion("Advanced Settings", open=False):
        with gr.Row():
            pbs_min_slider = gr.Slider(label="Min PBS Length", minimum=1, maximum=20, value=8, step=1)
            pbs_max_slider = gr.Slider(label="Max PBS Length", minimum=1, maximum=30, value=17, step=1)
        with gr.Row():
            rt_min_slider = gr.Slider(label="Min RT Overhang Length", minimum=1, maximum=20, value=1, step=1)
            rt_max_slider = gr.Slider(label="Max RT Overhang Length", minimum=1, maximum=40, value=31, step=1)

    run_button = gr.Button("Run Design", variant="primary")
    
    status_output = gr.Textbox(label="Status")
    
    with gr.Tab("pegRNA Results"):
        peg_output = gr.DataFrame(label="pegRNA Candidates")
        
    with gr.Tab("Nicking Guide Results"):
        nick_output = gr.DataFrame(label="Nicking Guide Candidates")

    run_button.click(
        fn=run_pridict,
        inputs=[sequence_input, editor_dropdown, pbs_min_slider, pbs_max_slider, rt_min_slider, rt_max_slider],
        outputs=[peg_output, nick_output, status_output]
    )

demo.launch()
