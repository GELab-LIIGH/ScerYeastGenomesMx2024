






LOG_FILE="FST_SACE_881100.err"

# Archivo de salida con los resultados organizados
OUTPUT_FILE="fst_summary_3.txt"

# Escribir encabezado en la tabla de salida
echo -e "Comparación\tMean_Fst\tWeighted_Fst" > "$OUTPUT_FILE"

# Variables para almacenar los valores
pop1=""
pop2=""
mean_fst=""
weighted_fst=""

# Leer el archivo línea por línea
while IFS= read -r line; do
    # Detectar nombres de poblaciones
    if [[ $line == *"--weir-fst-pop"* ]]; then
        if [[ -z "$pop1" ]]; then
            pop1=$(basename "$line" .txt | awk '{print $NF}')
        else
            pop2=$(basename "$line" .txt | awk '{print $NF}')
        fi
    fi

    # Detectar el Mean Fst
    if [[ $line == *"Weir and Cockerham mean Fst estimate"* ]]; then
        mean_fst=$(echo "$line" | awk '{print $NF}')
    fi

    # Detectar el Weighted Fst
    if [[ $line == *"Weir and Cockerham weighted Fst estimate"* ]]; then
        weighted_fst=$(echo "$line" | awk '{print $NF}')

        # Guardar el resultado cuando ambos valores han sido leídos
        if [[ -n "$pop1" && -n "$pop2" && -n "$mean_fst" && -n "$weighted_fst" ]]; then
            echo -e "${pop1}_vs_${pop2}\t${mean_fst}\t${weighted_fst}" >> "$OUTPUT_FILE"
        fi

        # Resetear variables para la siguiente comparación
        pop1=""
        pop2=""
        mean_fst=""
        weighted_fst=""
    fi
done < "$LOG_FILE"

echo "Tabla generada en $OUTPUT_FILE"
