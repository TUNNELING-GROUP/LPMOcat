import numpy as np
import plotly.graph_objects as go

def sec_struct_plot(Figure=None, sequence="", row=1, protein_id=None, three_letter_sequence=None, auth_numbers=None,
                    glyph_scale=1.0):
    """
    Elegant 1D visualization of protein secondary structure with shadows
    for both alpha helices and beta sheets.

    glyph_scale multiplies the heights of the helix waves, strand arrows and coil wiggle (drawn
    around y=0): 1 suits an autoranged y-axis, larger values a fixed one (see makePlots).
    """

    # Colors
    color_map = {
        'H': '#C62828',  # red
        'E': '#FBC02D',  # yellow
        'U': '#616161',  # grey
    }

    sec_structure_legens = {
        'H': 'Helix',
        'E': 'Sheet',
        'U': 'Undefined'
    }

    hoverList = []
    baseline = 0  # vertical track per protein
    pos = 0

    if Figure is None:
        Figure = go.Figure()

    legend_flags = {k: True for k in color_map.keys()}
    
    # Residue numbering per position. Prefer the structure's original (auth) residue
    # numbers so the SS plot matches the Sequence/Charge plots and the 3D viewer; fall
    # back to a 1-based count only if none were supplied.
    if auth_numbers is not None:
        residue_numbers = [auth_numbers[i] if i < len(auth_numbers) else None
                           for i in range(len(sequence))]
    else:
        residue_numbers = []
        residue_count = 1
        for i in range(len(sequence)):
            if sequence[i] == '-':
                residue_numbers.append(None)
            else:
                residue_numbers.append(residue_count)
                residue_count += 1

    # OPTIMIZATION: Pre-calculate gap positions to avoid repeated checks
    gap_positions = set(i for i, char in enumerate(sequence) if char == '-')
    
    while pos < len(sequence):
        char = sequence[pos]
        start = pos
        while pos < len(sequence) and sequence[pos] == char:
            pos += 1
        end = pos

        if char == "-":
            # Add one None per gap position to maintain alignment
            hoverList.extend([None] * (end - start))

        else:
            # OPTIMIZATION: Use set lookup instead of string comparison
            next_char_is_gap = (end in gap_positions) if end < len(sequence) else False
            next_char = sequence[end] if end < len(sequence) and not next_char_is_gap else None
            end_offset = 0.5 if next_char_is_gap else 1.0
            
            # OPTIMIZATION: Use set lookup instead of string comparison
            prev_char_was_gap = (start - 1 in gap_positions) if start > 0 else False
            prev_char = sequence[start - 1] if start > 0 and not prev_char_was_gap else None
            start_offset = -0.5 if prev_char_was_gap else 0.0

            showlegend = legend_flags[char]
            legend_flags[char] = False

            # Calculate start and end residue numbers for this structure element (excluding gaps)
            start_residue = residue_numbers[start] if residue_numbers[start] is not None else None
            end_residue = residue_numbers[end-1] if residue_numbers[end-1] is not None else None

            # OPTIMIZATION: Pre-build hover text template to avoid repeated string operations
            struct_name = sec_structure_legens[char]
            if start_residue is not None and end_residue is not None:
                range_text = f" ({start_residue}-{end_residue})"
            else:
                range_text = ""
            
            # Add individual hover points for each residue in this structure element
            for i in range(start, end):
                residue_num = residue_numbers[i]
                if three_letter_sequence and i < len(three_letter_sequence):
                    residue_name = three_letter_sequence[i]
                    if residue_num is not None:
                        hover_text = f"{residue_name}{residue_num}<br>Secondary structure: {struct_name}{range_text}"
                    else:
                        hover_text = f"{residue_name}<br>Secondary structure: {struct_name}"
                else:
                    if residue_num is not None:
                        hover_text = f"Position {residue_num}<br>Secondary structure: {struct_name}"
                    else:
                        hover_text = None

                hoverList.append(hover_text)

                # Add invisible hover point for each residue
                # Figure.add_trace(
                #     go.Scatter(
                #         x=[i + 1.5],  # Center of residue position
                #         y=[baseline],
                #         mode='markers',
                #         marker=dict(size=0.1, color='rgba(0,0,0,0)'),  # Invisible
                #         showlegend=False,
                #         hovertext=hover_text,
                #         hoverinfo='text'
                #     )
                # )

            # --- Alpha Helix ---
            if char == 'H':
                segment_length = end - start
                
                # Unified rule: all segments fill cursor box from start+0.5 to end+0.5
                if segment_length == 1:
                    # Single residue: fill cursor box
                    exact_position = start + 1
                    num_points = 10
                    helix_x = np.linspace(exact_position - 0.5, exact_position + 0.5, num_points)
                    res_per_wave = 3
                    waves = max(1, segment_length / res_per_wave)
                    phase = (helix_x - (exact_position - 0.5)) * 2 * np.pi * waves / segment_length
                    helix_y = 0.10 * glyph_scale * np.sin(phase) + baseline
                    
                    Figure.add_trace(
                        go.Scatter(
                            x=helix_x, y=helix_y,
                            mode='lines',
                            line=dict(color=color_map[char], width=4, shape='spline'),
                            name='Helix',
                            showlegend=showlegend,
                            hoverinfo='skip'
                        )
                    )
                else:
                    # Multi-residue: fill cursor box from start+0.5 to end+0.5
                    if segment_length < 10:
                        num_points = max(15, segment_length * 3)
                    elif segment_length < 20:
                        num_points = segment_length * 4
                    else:
                        num_points = max(40, segment_length * 4)
                    
                    helix_start = start + 0.5
                    helix_end = end + 0.5
                    helix_x = np.linspace(helix_start, helix_end, num_points)
                    res_per_wave = 3
                    waves = max(1, segment_length / res_per_wave)
                    phase = (helix_x - helix_start) * 2 * np.pi * waves / segment_length
                    helix_y = 0.10 * glyph_scale * np.sin(phase) + baseline

                    Figure.add_trace(
                        go.Scatter(
                            x=helix_x, y=helix_y,
                            mode='lines',
                            line=dict(color=color_map[char], width=4, shape='spline'),
                            name='Helix',
                            showlegend=showlegend,
                            hoverinfo='skip'
                        )
                    )

            # --- Beta Sheet ---
            elif char == 'E':
                segment_length = end - start
                body_height = 0.015 * glyph_scale
                arrow_height = 0.05 * glyph_scale
                
                # Unified rule: all segments fill cursor box from start+0.5 to end+0.5
                x0 = start + 0.5
                x1 = end + 0.5
                
                if segment_length == 1:
                    # Single residue: rectangle only, no arrowhead
                    Figure.add_shape(
                        type="rect",
                        x0=x0, y0=baseline-body_height,
                        x1=x1, y1=baseline+body_height,
                        line=dict(width=0),
                        fillcolor=color_map[char],
                        xref='x', yref='y'
                    )
                else:
                    # Multi-residue: rectangle body + arrowhead
                    arrow_width = min(0.8, segment_length)
                    # Main rectangle body - fill cursor box minus arrow
                    rect_x1 = end + 0.5 - arrow_width
                    
                    Figure.add_shape(
                        type="rect",
                        x0=x0, y0=baseline-body_height,
                        x1=rect_x1, y1=baseline+body_height,
                        line=dict(width=0),
                        fillcolor=color_map[char],
                        xref='x', yref='y'
                    )
                    # Main arrowhead - end at end+0.5
                    arrow_x = end + 0.5
                    Figure.add_shape(
                        type="path",
                        path=f"M {arrow_x-arrow_width},{baseline-arrow_height} L {arrow_x},{baseline} L {arrow_x-arrow_width},{baseline+arrow_height} Z",
                        fillcolor=color_map[char],
                        line=dict(width=0),
                        xref='x', yref='y'
                    )
                
                # Invisible scatter for legend only
                Figure.add_trace(
                    go.Scatter(
                        x=[(start+end)/2 + 1], y=[baseline],
                        mode='markers',
                        marker=dict(size=0.1, color=color_map[char]),
                        name='Sheet',
                        showlegend=showlegend,
                        hoverinfo='skip'
                    )
                )

            # --- Undefined ---
            elif char == 'U':
                segment_length = end - start
                
                # Unified rule: all segments fill cursor box from start+0.5 to end+0.5
                if segment_length < 3:
                    num_points = segment_length * 2
                elif segment_length < 10:
                    num_points = segment_length * 3
                else:
                    num_points = max(25, segment_length * 4)
                
                coil_start = start + 0.5
                coil_end = end + 0.5
                coil_x = np.linspace(coil_start, coil_end, num_points)
                phase_base = coil_start

                # --- CONSTANT FREQUENCY ---
                res_per_wave = 3   # same unit as helix
                waves = segment_length / res_per_wave
                amplitude = 0.01 * glyph_scale   # vertical amplitude for loop
                coil_y = amplitude * np.sin((coil_x - phase_base) * 2 * np.pi * waves / segment_length) + baseline

                Figure.add_trace(
                    go.Scatter(
                        x=coil_x, y=coil_y,
                        mode='lines',
                        line=dict(color=color_map[char], width=3, shape='spline'),
                        name='Undefined',
                        showlegend=showlegend,
                        hoverinfo='skip'
                    )
                )


    # Ensure full sequence is visible
    Figure.add_trace(
        go.Scatter(
            x=[1, len(sequence) + 1], y=[baseline, baseline],  # Shift by +1
            mode='markers',
            marker=dict(size=0.1, color='white'),
            showlegend=False,
            hoverinfo='skip'
        )
    )

    Figure.update_layout(
    plot_bgcolor='white',
    xaxis=dict(
        showgrid=False, zeroline=False, showline=False, ticks='',
        range=[1, len(sequence) + 1]   # <-- forces x axis to start at 1
    ),
    yaxis=dict(showgrid=False, zeroline=False, showline=False, showticklabels=False),
    margin=dict(l=40, r=20, t=20, b=20)
    )


    # adding trace that makes plots aligned with proper hover handling
    # Replace None values with empty string to prevent "null" from showing
    hover_data = [h if h is not None else '' for h in hoverList]
    
    # Convert hex colors to rgba for hover backgrounds
    def hex_to_rgba(hex_color, alpha=0.7):
        """Convert hex color to rgba format"""
        hex_color = hex_color.lstrip('#')
        r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
        return f'rgba({r},{g},{b},{alpha})'
    
    # Create list of hover background colors based on secondary structure
    hover_colors = []
    for char in sequence:
        if char == '-':
            # Use light gray for gaps
            hover_colors.append('rgba(176,190,197,0.7)')
        elif char in color_map:
            hover_colors.append(hex_to_rgba(color_map[char]))
        else:
            # Default to light gray for unknown characters
            hover_colors.append('rgba(176,190,197,0.7)')
    
    # Filter out gap positions for hover (only show hovers for actual residues)
    x_no_gaps = []
    y_no_gaps = []
    customdata_no_gaps = []
    hover_colors_no_gaps = []
    for idx, hd in enumerate(hover_data):
        if hd != '':  # Skip gaps (empty strings)
            x_no_gaps.append(idx + 1)
            y_no_gaps.append(baseline)
            customdata_no_gaps.append(hd)
            hover_colors_no_gaps.append(hover_colors[idx])
    
    Figure.add_trace(go.Scatter(
        x=x_no_gaps,
        y=y_no_gaps,
        line_width=0, 
        customdata=customdata_no_gaps, 
        hovertemplate='%{customdata}<extra></extra>',
        hoverinfo='text',
        hoverlabel=dict(bgcolor=hover_colors_no_gaps),
    ))

    return Figure
